# Copyright (C) 2004-2006 Paul B. Horton
# Copyright 2026 Google LLC
# Licensed under the GNU General Public License v3 (GPLv3).

package WoLFPSort::KNN;

use strict;
use warnings;
use Carp         qw(croak);
use Scalar::Util qw(blessed);
use WoLFPSort;
use WoLFPSort::ModelLoader;
use WoLFPSort::UniformScaler;

our $VERSION = '0.200';

my $_PDL_LOADED = 0;

sub _ensure_pdl {
    return if $_PDL_LOADED;
    require PDL::Lite;
    require PDL::Core;
    require PDL::Ufunc;
    $_PDL_LOADED = 1;
    return;
}

sub new {
    my ( $class, %args ) = @_;
    for my $req (qw(weight_file training_file utility_file)) {
        croak "Missing required argument '$req'" unless defined $args{$req};
    }

    my $use_xs = WoLFPSort::has_xs();
    if ($use_xs) {
        my $xs
            = WoLFPSort::XS->new( $args{weight_file}, $args{training_file}, $args{utility_file}, );
        return bless {
            use_xs  => 1,
            use_pdl => 0,
            xs      => $xs,
            k       => $xs->k(),
        }, $class;
    }

    my $use_pdl = WoLFPSort::has_pdl();
    _ensure_pdl() if $use_pdl;

    my $w = WoLFPSort::ModelLoader::load_weights( $args{weight_file} );
    my $u = WoLFPSort::ModelLoader::load_utility( $args{utility_file} );
    my $t = WoLFPSort::ModelLoader::load_training( $args{training_file}, $w->{feature_names},
        $u->{class_index}, pack_doubles => $use_pdl, );

    my ( $scaler, $scaled_X, $class_ids );
    if ($use_pdl) {
        my $raw_X = PDL::Core::zeroes( $t->{n_feat}, $t->{n_rows} );
        ${ $raw_X->get_dataref } = $t->{packed_buf};
        $raw_X->upd_data;
        $scaler    = WoLFPSort::UniformScaler->new($raw_X);
        $scaled_X  = $scaler->transform($raw_X);
        $class_ids = PDL::Core::long( $t->{class_ids} );
    }
    else {
        $scaler    = WoLFPSort::UniformScaler->new( $t->{rows} );
        $scaled_X  = $scaler->transform( $t->{rows} );
        $class_ids = $t->{class_ids};
    }

    return bless {
        use_pdl       => $use_pdl,
        feature_names => $w->{feature_names},
        w1            => $use_pdl ? PDL::Core::double( $w->{w1} ) : $w->{w1},
        w2            => $use_pdl ? PDL::Core::double( $w->{w2} ) : $w->{w2},
        k             => $w->{k},
        class_names   => $u->{class_names},
        class_index   => $u->{class_index},
        utility       => $use_pdl ? PDL::Core::double( $u->{matrix} ) : $u->{matrix},
        ids           => $t->{ids},
        classes       => $t->{classes},
        class_ids     => $class_ids,
        scaler        => $scaler,
        scaled_X      => $scaled_X,
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
    return $self->{xs}->n_samples() if $self->{use_xs};
    return scalar @{ $self->{ids} };
}

sub _predict_pp {
    my ( $self, $query_vec, $query_id, $ignore_self ) = @_;
    my $q_scaled = $self->{scaler}->transform($query_vec);
    my $X        = $self->{scaled_X};
    my $w1       = $self->{w1};
    my $w2       = $self->{w2};
    my $n        = $self->n_samples;
    my $n_feat   = scalar @{$w1};

    my @dists;
    for my $i ( 0 .. $n - 1 ) {
        my $row = $X->[$i];
        my ( $l1, $l2sq ) = ( 0.0, 0.0 );
        for my $f ( 0 .. $n_feat - 1 ) {
            my $d = $row->[$f] - $q_scaled->[$f];
            $l1   += $w1->[$f] * abs $d  if $w1->[$f];
            $l2sq += $w2->[$f] * $d * $d if $w2->[$f];
        }
        push @dists, $l1 + sqrt $l2sq;
    }

    my @order = sort { $dists[$a] <=> $dists[$b] } ( 0 .. $n - 1 );
    my $k     = $self->{k};
    $k = $n if $k > $n;

    my @neighbor_indices;
    my $ids = $self->{ids};
    for my $idx (@order) {
        if ( $ignore_self && defined $query_id && $ids->[$idx] eq $query_id ) {
            next;
        }
        push @neighbor_indices, $idx;
        last if @neighbor_indices == $k;
    }

    my $n_cls  = scalar @{ $self->{class_names} };
    my @counts = (0.0) x $n_cls;
    for my $idx (@neighbor_indices) {
        $counts[ $self->{class_ids}[$idx] ] += 1.0;
    }

    my $util = $self->{utility};
    my @scores;
    for my $c ( 0 .. $n_cls - 1 ) {
        my $s = 0.0;
        for my $j ( 0 .. $n_cls - 1 ) {
            $s += $util->[$c][$j] * $counts[$j];
        }
        push @scores, $s;
    }

    my @nb_dists = map { $dists[$_] } @neighbor_indices;
    return $self->_format_result( \@scores, \@neighbor_indices, \@nb_dists );
}

sub predict_vector {
    my ( $self, $query_vec, %opts ) = @_;
    my $query_id    = $opts{id};
    my $ignore_self = $opts{ignore_self} ? 1 : 0;

    if ( $self->{use_xs} ) {
        my $arr
            = blessed($query_vec) && $query_vec->isa('PDL')
            ? [ $query_vec->list ]
            : $query_vec;
        return $self->{xs}->predict_vector( $arr, $query_id, $ignore_self );
    }

    if ( !$self->{use_pdl} ) {
        return $self->_predict_pp( $query_vec, $query_id, $ignore_self );
    }

    my $q_pdl
        = blessed($query_vec) && $query_vec->isa('PDL')
        ? $query_vec
        : PDL::Core::double($query_vec);

    my $q_scaled = $self->{scaler}->transform($q_pdl);
    my $diff     = $self->{scaled_X} - $q_scaled;

    # Horton et al. (APBC06), Section 2.3.2, Equation 1:
    # d_W(j, k) = sum_i w1_i * |F_ji - F_ki| + sqrt(sum_i w2_i * (F_ji - F_ki)^2)
    my $l1   = PDL::Ufunc::sumover( $self->{w1} * abs $diff );
    my $l2   = sqrt PDL::Ufunc::sumover( $self->{w2} * ( $diff * $diff ) );
    my $dist = $l1 + $l2;

    my $order = PDL::Ufunc::qsorti($dist);
    my $k     = $self->{k};
    my $n     = $self->n_samples;
    $k = $n if $k > $n;

    my @neighbor_indices;
    my $ids = $self->{ids};
    for my $pos ( 0 .. $n - 1 ) {
        my $idx = $order->at($pos);
        if ( $ignore_self && defined $query_id && $ids->[$idx] eq $query_id ) {
            next;
        }
        push @neighbor_indices, $idx;
        last if @neighbor_indices == $k;
    }

    my $n_cls  = scalar @{ $self->{class_names} };
    my $counts = PDL::Core::zeroes($n_cls);
    for my $idx (@neighbor_indices) {
        my $c = $self->{class_ids}->at($idx);
        $counts->set( $c, $counts->at($c) + 1.0 );
    }

    # Horton et al. (APBC06), Section 2.3.3, Table 3:
    # Expected utility score for each candidate site using the symmetric utility matrix.
    my @scores   = PDL::Ufunc::sumover( $self->{utility} * $counts )->list;
    my @nb_dists = map { $dist->at($_) } @neighbor_indices;
    return $self->_format_result( \@scores, \@neighbor_indices, \@nb_dists );
}

sub _format_result {
    my ( $self, $scores, $neighbor_indices, $nb_dists ) = @_;
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

Weighted hybrid L1/L2 I<k>-nearest neighbors classifier replacing the 32-bit
C<wolfPredict> binary. Uses L<PDL> when available and falls back to pure Perl
otherwise.

=head1 TYPES

When L<PDL> is available, feature matrices, weights, and utility matrices are
C<double> ndarrays and class indices are C<long>. Complex types are not
supported.

=cut
