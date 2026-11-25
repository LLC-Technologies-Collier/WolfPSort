# Copyright (C) 2004-2006 Paul B. Horton
# Copyright 2026 Google LLC
# Licensed under the GNU General Public License v3 (GPLv3).

package WoLFPSort::PSORTII::KNN;

use strict;
use warnings;
use Carp qw(croak);
use WoLFPSort::Backend;

our $VERSION = '0.200';

sub new {
    my ( $class, %args ) = @_;
    for my $req (qw(matrix names locations)) {
        croak "Missing required argument '$req'" unless defined $args{$req};
    }

    my $backend = WoLFPSort::Backend->current();
    my $mat     = $backend->psortii_prepare( $args{matrix} );

    my %class_freq;
    for my $loc ( @{ $args{locations} } ) {
        $class_freq{$loc}++;
    }

    return bless {
        backend    => $backend,
        matrix     => $mat,
        names      => $args{names},
        locations  => $args{locations},
        class_freq => \%class_freq,
    }, $class;
}

sub predict {
    my ( $self, $query, $k ) = @_;
    $k ||= 25;

    my $d    = $self->{backend}->psortii_distances( $self->{matrix}, $query );
    my $locs = $self->{locations};
    my $cf   = $self->{class_freq};
    my $n    = scalar @{ $self->{names} };
    $k = $n if $k > $n;

    my @order = sort { $d->[$a] <=> $d->[$b] || $cf->{ $locs->[$a] } <=> $cf->{ $locs->[$b] } }
        ( 0 .. $n - 1 );

    my %freq;
    for my $i ( 0 .. $k - 1 ) {
        my $loc = $locs->[ $order[$i] ];
        $freq{$loc}++;
    }

    my @sorted_classes = sort { $freq{$b} <=> $freq{$a} || $cf->{$b} <=> $cf->{$a} } keys %freq;

    return {
        prediction => $sorted_classes[0],
        classes    => [
            map { +{ class => $_, count => $freq{$_}, percent => ( $freq{$_} / $k ) * 100.0, } }
                @sorted_classes
        ],
    };
}

1;

__END__

=head1 NAME

WoLFPSort::PSORTII::KNN - Legacy PSORT II Euclidean kNN classifier

=head1 SYNOPSIS

    use WoLFPSort::PSORTII::KNN;

    my $knn = WoLFPSort::PSORTII::KNN->new(
        matrix    => \@rows,
        names     => \@seq_ids,
        locations => \@loc_labels,
    );
    my $res = $knn->predict(\@query_features, 25);

=cut
