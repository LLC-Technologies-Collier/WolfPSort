# Copyright (C) 2004-2006 Paul B. Horton
# Copyright 2026 Google LLC
# Licensed under the GNU General Public License v3 (GPLv3).

package WoLFPSort::Backend::PurePerl;

use strict;
use warnings;
use Carp qw(croak);

our $VERSION = '0.200';

sub is_pdl {
    return 0;
}

sub build_scaler_tables {
    my ( $class, $rows ) = @_;
    croak 'Expected a 2D PDL ndarray (F x N) or arrayref of feature rows'
        unless ref($rows) eq 'ARRAY' && @{$rows} && ref( $rows->[0] ) eq 'ARRAY';
    my $n_samples  = scalar @{$rows};
    my $n_features = scalar @{ $rows->[0] };
    croak 'Training matrix must have at least 1 feature and 1 sample'
        if $n_features < 1 || $n_samples < 1;

    my @tables;
    for my $f ( 0 .. $n_features - 1 ) {
        my @col = sort { $a <=> $b } map { $_->[$f] + 0 } @{$rows};
        my ( @vals, @lo, @ranks );
        my $i = 0;
        while ( $i < $n_samples ) {
            my $v = $col[$i];
            my $j = $i + 1;
            $j++ while $j < $n_samples && $col[$j] == $v;
            push @vals, $v;
            push @lo,   $i;
            push @ranks, ( 0.5 * ( $i + $j - 1.0 ) / $n_samples ) * 100.0;
            $i = $j;
        }
        push @tables,
            {
            vals  => \@vals,
            lo    => \@lo,
            ranks => \@ranks,
            n_u   => scalar @vals,
            };
    }
    return ( \@tables, $n_features, $n_samples );
}

sub _bsearch_left {
    my ( $x, $vals, $n_u ) = @_;
    my ( $lo, $hi ) = ( 0, $n_u );
    while ( $lo < $hi ) {
        my $mid = int( ( $lo + $hi ) / 2 );
        if ( $vals->[$mid] < $x ) {
            $lo = $mid + 1;
        }
        else {
            $hi = $mid;
        }
    }
    return $lo;
}

# Horton et al. (APBC06), Section 2.3.1 (Candidate Features), p. 3, para. 1,
# sentences 3-4: linear interpolation (and extrapolation outside [min, max]).
sub _interp {
    my ( $x, $vals, $ranks, $n_u, $idx ) = @_;
    return $ranks->[0] if $n_u <= 1;
    my $i0 = $idx - 1;
    $i0 = 0        if $i0 < 0;
    $i0 = $n_u - 2 if $i0 > $n_u - 2;
    my $i1 = $i0 + 1;
    my $dv = $vals->[$i1] - $vals->[$i0];
    return $ranks->[$i0] if $dv == 0;
    return $ranks->[$i0] + ( ( $x - $vals->[$i0] ) / $dv ) * ( $ranks->[$i1] - $ranks->[$i0] );
}

sub _transform_row {
    my ( $tables, $n_features, $row ) = @_;
    croak 'Feature count mismatch: got ' . scalar( @{$row} ) . ", expected $n_features"
        if @{$row} != $n_features;

    my @out;
    for my $f ( 0 .. $n_features - 1 ) {
        my $tbl = $tables->[$f];
        my $x   = $row->[$f] + 0;
        my $n_u = $tbl->{n_u};
        my $idx = _bsearch_left( $x, $tbl->{vals}, $n_u );
        if ( $idx < $n_u && $tbl->{vals}[$idx] == $x ) {
            push @out, $tbl->{ranks}[$idx];
        }
        else {
            push @out, _interp( $x, $tbl->{vals}, $tbl->{ranks}, $n_u, $idx );
        }
    }
    return \@out;
}

sub transform_scaler {
    my ( $class, $tables, $n_features, $X ) = @_;
    croak 'Expected an arrayref of features or arrayref of rows'
        unless ref($X) eq 'ARRAY';
    if ( @{$X} && ref( $X->[0] ) eq 'ARRAY' ) {
        return [ map { _transform_row( $tables, $n_features, $_ ) } @{$X} ];
    }
    return _transform_row( $tables, $n_features, $X );
}

sub prepare_knn_model {
    my ( $class, $w, $u, $t ) = @_;
    require WoLFPSort::UniformScaler;
    my $scaler = WoLFPSort::UniformScaler->new( $t->{rows} );
    return {
        w1        => $w->{w1},
        w2        => $w->{w2},
        utility   => $u->{matrix},
        class_ids => $t->{class_ids},
        scaler    => $scaler,
        scaled_X  => $scaler->transform( $t->{rows} ),
    };
}

sub predict_knn {
    my ( $class, $args ) = @_;
    my ( $state, $query_vec, $ids, $k, $n_cls, $query_id, $ignore_self )
        = @{$args}{qw(state query_vec ids k n_cls query_id ignore_self)};
    my $q_scaled = $state->{scaler}->transform($query_vec);
    my $X        = $state->{scaled_X};
    my $w1       = $state->{w1};
    my $w2       = $state->{w2};
    my $n        = scalar @{$ids};
    my $n_feat   = scalar @{$w1};

    # Horton et al. (APBC06), Section 2.3.2, p. 4, para. 1, Equation 1
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
    $k = $n if $k > $n;

    my @neighbor_indices;
    for my $idx (@order) {
        if ( $ignore_self && defined $query_id && $ids->[$idx] eq $query_id ) {
            next;
        }
        push @neighbor_indices, $idx;
        last if @neighbor_indices == $k;
    }

    my @counts = (0.0) x $n_cls;
    for my $idx (@neighbor_indices) {
        $counts[ $state->{class_ids}[$idx] ] += 1.0;
    }

    # Horton et al. (APBC06), Section 2.3.3, p. 4, para. 1 & Table 3
    my $util = $state->{utility};
    my @scores;
    for my $c ( 0 .. $n_cls - 1 ) {
        my $s = 0.0;
        for my $j ( 0 .. $n_cls - 1 ) {
            $s += $util->[$c][$j] * $counts[$j];
        }
        push @scores, $s;
    }

    my @nb_dists = map { $dists[$_] } @neighbor_indices;
    return ( \@scores, \@neighbor_indices, \@nb_dists );
}

sub vector_sum {
    my ( $class, $vals ) = @_;
    my $sum = 0.0;
    for my $v ( @{$vals} ) {
        $sum += $v;
    }
    return $sum;
}

sub weighted_average {
    my ( $class, $weights, $vals ) = @_;
    my ( $sum, $wsum ) = ( 0.0, 0.0 );
    for my $i ( 0 .. $#{$vals} ) {
        $sum  += $weights->[$i] * $vals->[$i];
        $wsum += $weights->[$i];
    }
    return $sum / $wsum;
}

sub max_window_sum {
    my ( $class, $vals, $window_len, $n_win ) = @_;
    my $cur = 0.0;
    for my $i ( 0 .. $window_len - 1 ) {
        $cur += $vals->[$i];
    }
    my $max = $cur;
    for my $w ( 1 .. $n_win - 1 ) {
        $cur -= $vals->[ $w - 1 ];
        $cur += $vals->[ $w + $window_len - 1 ];
        $max = $cur if $cur > $max;
    }
    return $max;
}

sub encode_indices {
    my ( $class, $idx_ref ) = @_;
    return $idx_ref;
}

sub psortii_prepare {
    my ( $class, $matrix ) = @_;
    croak 'Expected arrayref of training rows' unless ref($matrix) eq 'ARRAY';
    return $matrix;
}

sub psortii_distances {
    my ( $class, $mat, $query ) = @_;
    my @dists;
    for my $row ( @{$mat} ) {
        my $dist = 0.0;
        for my $f ( 0 .. $#{$row} ) {
            my $delta = $row->[$f] - $query->[$f];
            $dist += $delta * $delta;
        }
        push @dists, $dist;
    }
    return \@dists;
}

1;

__END__

=head1 NAME

WoLFPSort::Backend::PurePerl - Pure-Perl numerical backend for WoLFPSort

=cut
