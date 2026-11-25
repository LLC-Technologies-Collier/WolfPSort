# Copyright (C) 2004-2006 Paul B. Horton
# Copyright 2026 Google LLC
# Licensed under the GNU General Public License v3 (GPLv3).

package WoLFPSort::KNN;

use strict;
use warnings;
use Carp         qw(croak);
use Scalar::Util qw(blessed);
use WoLFPSort;
use WoLFPSort::Backend;
use WoLFPSort::ModelLoader;

our $VERSION = '0.200';

sub new {
    my ( $class, %args ) = @_;
    for my $req (qw(weight_file training_file utility_file)) {
        croak "Missing required argument '$req'" unless defined $args{$req};
    }

    if ( WoLFPSort::has_xs() ) {
        my $xs
            = WoLFPSort::XS->new( $args{weight_file}, $args{training_file}, $args{utility_file}, );
        return bless {
            xs => $xs,
            k  => $xs->k(),
        }, $class;
    }

    my $backend = WoLFPSort::Backend->current();
    my $w       = WoLFPSort::ModelLoader::load_weights( $args{weight_file} );
    my $u       = WoLFPSort::ModelLoader::load_utility( $args{utility_file} );
    my $t       = WoLFPSort::ModelLoader::load_training( $args{training_file}, $w->{feature_names},
        $u->{class_index}, pack_doubles => $backend->is_pdl(), );

    my $state = $backend->prepare_knn_model( $w, $u, $t );

    return bless {
        backend       => $backend,
        state         => $state,
        feature_names => $w->{feature_names},
        k             => $w->{k},
        class_names   => $u->{class_names},
        ids           => $t->{ids},
        classes       => $t->{classes},
    }, $class;
}

sub k {
    my ($self) = @_;
    return $self->{k};
}

sub feature_names {
    my ($self) = @_;
    return $self->{feature_names};
}

sub class_names {
    my ($self) = @_;
    return $self->{class_names};
}

sub n_samples {
    my ($self) = @_;
    return $self->{xs}->n_samples() if $self->{xs};
    return scalar @{ $self->{ids} };
}

sub predict_vector {
    my ( $self, $query_vec, %opts ) = @_;
    my $query_id    = $opts{id};
    my $ignore_self = $opts{ignore_self} ? 1 : 0;

    if ( $self->{xs} ) {
        my $arr
            = blessed($query_vec) && $query_vec->isa('PDL')
            ? [ $query_vec->list ]
            : $query_vec;
        return $self->{xs}->predict_vector( $arr, $query_id, $ignore_self );
    }

    my ( $scores, $neighbor_indices, $nb_dists ) = $self->{backend}->predict_knn(
        {   state       => $self->{state},
            query_vec   => $query_vec,
            ids         => $self->{ids},
            k           => $self->{k},
            n_cls       => scalar( @{ $self->{class_names} } ),
            query_id    => $query_id,
            ignore_self => $ignore_self,
        }
    );

    my $n_cls  = scalar @{ $self->{class_names} };
    my @ranked = sort { $scores->[$b] <=> $scores->[$a] || $a <=> $b } ( 0 .. $n_cls - 1 );

    my @all_scores
        = map { +{ class => $self->{class_names}[$_], score => $scores->[$_], } } @ranked;

    my @neighbors;
    for my $i ( 0 .. $#{$neighbor_indices} ) {
        my $idx = $neighbor_indices->[$i];
        push @neighbors,
            {
            id         => $self->{ids}[$idx],
            class      => $self->{classes}[$idx],
            similarity => -$nb_dists->[$i],
            };
    }

    return {
        scores    => \@all_scores,
        neighbors => \@neighbors,
    };
}

sub filter_scores {
    my ( $self, $scores, $print_all ) = @_;
    require WoLFPSort::CLI::Predict;
    return WoLFPSort::CLI::Predict::filter_scores( $scores, $print_all );
}

sub format_number {
    my ($val) = @_;
    require WoLFPSort::CLI::Predict;
    return WoLFPSort::CLI::Predict::format_number($val);
}

1;

__END__

=head1 NAME

WoLFPSort::KNN - Weighted hybrid L1/L2 kNN classifier for WoLF PSORT

=head1 SYNOPSIS

    use WoLFPSort::KNN;

    my $knn = WoLFPSort::KNN->new(
        weight_file   => 'data/fungi.wolfw',
        training_file => 'data/fungi.wolff',
        utility_file  => 'data/fungi.wolfu',
    );

    my $res = $knn->predict_vector(\@features, id => 'RCC1_YEAST');

=head1 DESCRIPTION

Weighted hybrid L1/L2 I<k>-nearest neighbors classifier (Horton et al., APBC06,
Sections 2.3.2-2.3.3). Dispatches to L<WoLFPSort::XS> when compiled, or to
L<WoLFPSort::Backend> (PDL or pure Perl).

=cut
