# Copyright (C) 2004-2006 Paul B. Horton
# Copyright 2026 Google LLC
# Licensed under the GNU General Public License v3 (GPLv3).

package WoLFPSort::PSORTII::KNN;

use strict;
use warnings;
use Carp         qw(croak);
use Scalar::Util qw(blessed);
use WoLFPSort;

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
    for my $req (qw(matrix names locations)) {
        croak "Missing required argument '$req'" unless defined $args{$req};
    }

    my $use_pdl = WoLFPSort::has_pdl();
    my $mat;
    if ($use_pdl) {
        _ensure_pdl();
        $mat
            = blessed( $args{matrix} ) && $args{matrix}->isa('PDL')
            ? $args{matrix}
            : PDL::Core::double( $args{matrix} );
        croak 'Expected 2D training matrix (F x N)' unless $mat->ndims == 2;
    }
    else {
        croak 'Expected arrayref of training rows'
            unless ref( $args{matrix} ) eq 'ARRAY';
        $mat = $args{matrix};
    }

    my %class_freq;
    for my $loc ( @{ $args{locations} } ) {
        $class_freq{$loc}++;
    }

    return bless {
        use_pdl    => $use_pdl,
        matrix     => $mat,
        names      => $args{names},
        locations  => $args{locations},
        class_freq => \%class_freq,
    }, $class;
}

sub predict {
    my ( $self, $query, $k ) = @_;
    $k ||= 25;

    my @d;
    if ( $self->{use_pdl} ) {
        my $q
            = blessed($query) && $query->isa('PDL')
            ? $query
            : PDL::Core::double($query);
        my $diff = $self->{matrix} - $q;
        @d = PDL::Ufunc::sumover( $diff * $diff )->list;
    }
    else {
        for my $row ( @{ $self->{matrix} } ) {
            my $dist = 0.0;
            for my $f ( 0 .. $#{$row} ) {
                my $delta = $row->[$f] - $query->[$f];
                $dist += $delta * $delta;
            }
            push @d, $dist;
        }
    }

    my $locs = $self->{locations};
    my $cf   = $self->{class_freq};
    my $n    = scalar @{ $self->{names} };
    $k = $n if $k > $n;

    my @order = sort { $d[$a] <=> $d[$b] || $cf->{ $locs->[$a] } <=> $cf->{ $locs->[$b] } }
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

=head1 TYPES

When L<PDL> is available, feature matrices and query vectors use C<double>
ndarrays. Complex types are not supported.

=cut
