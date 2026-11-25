# Copyright (C) 2004-2006 Paul B. Horton
# Copyright 2026 Google LLC
# Licensed under the GNU General Public License v3 (GPLv3).

package WoLFPSort::Backend::PDL;

use strict;
use warnings;
use Carp         qw(croak);
use Scalar::Util qw(blessed);
use PDL::Lite;
use PDL::Core      qw(double long zeroes);
use PDL::Ufunc     qw(cumusumover max qsorti sumover);
use PDL::Slices    qw(rle);
use PDL::Primitive qw(vsearch which);

our $VERSION = '0.200';

sub is_pdl {
    return 1;
}

sub build_scaler_tables {
    my ( $class, $matrix ) = @_;
    croak 'Expected a 2D PDL ndarray (F x N)'
        unless blessed($matrix) && $matrix->isa('PDL') && $matrix->ndims == 2;
    my ( $n_features, $n_samples ) = $matrix->dims;
    croak 'Training matrix must have at least 1 feature and 1 sample'
        if $n_features < 1 || $n_samples < 1;

    my $sorted = $matrix->transpose->qsort->sever;
    my @tables;
    for my $f ( 0 .. $n_features - 1 ) {
        my $s_col = $sorted->slice(":,($f)");
        my ( $counts, $vals ) = rle($s_col);
        my $nz = which( $counts > 0 );
        my $c  = $counts->index($nz)->double;
        my $v  = $vals->index($nz)->sever;
        my $hi = cumusumover($c);
        my $lo = $hi - $c;
        my $r  = ( 0.5 * ( $lo + $hi - 1.0 ) / $n_samples ) * 100.0;
        push @tables,
            {
            vals  => $v,
            lo    => $lo,
            ranks => $r,
            n_u   => $v->nelem,
            };
    }
    return ( \@tables, $n_features, $n_samples );
}

sub transform_scaler {
    my ( $class, $tables, $n_features, $X ) = @_;
    croak 'Expected a 1D or 2D PDL ndarray'
        unless blessed($X) && $X->isa('PDL');

    my $is_1d = ( $X->ndims == 1 );
    my $mat   = $is_1d ? $X->dummy(1) : $X;
    croak 'Expected 1D (F) or 2D (F x M) PDL ndarray'
        unless $mat->ndims == 2;

    my ( $n_feat, $n_obs ) = $mat->dims;
    croak "Feature count mismatch: got $n_feat, expected $n_features"
        if $n_feat != $n_features;

    my $mat_t = $mat->transpose->sever;
    my $out_t = zeroes( $n_obs, $n_feat );

    for my $f ( 0 .. $n_feat - 1 ) {
        my $tbl   = $tables->[$f];
        my $v     = $tbl->{vals};
        my $ranks = $tbl->{ranks};
        my $n_u   = $tbl->{n_u};

        my $x_col = $mat_t->slice(":,($f)");
        if ( $n_u <= 1 ) {
            $out_t->slice(":,($f)") .= $ranks->at(0);
            next;
        }

        # Horton et al. (APBC06), Section 2.3.1, p. 3, para. 1, sentences 3-4:
        # exact mid-rank for seen values, linear interpolation/extrapolation for unseen.
        my $idx       = vsearch( $x_col, $v, { mode => 'insert_leftmost' } );
        my $in_bounds = ( $idx < $n_u );
        my $safe_idx  = $idx->clip( 0, $n_u - 1 );
        my $is_match  = $in_bounds & ( $v->index($safe_idx) == $x_col );

        my $i0     = ( $idx - 1 )->clip( 0, $n_u - 2 );
        my $i1     = $i0 + 1;
        my $v0     = $v->index($i0);
        my $v1     = $v->index($i1);
        my $r0     = $ranks->index($i0);
        my $r1     = $ranks->index($i1);
        my $interp = $r0 + ( ( $x_col - $v0 ) / ( $v1 - $v0 ) ) * ( $r1 - $r0 );

        $out_t->slice(":,($f)") .= $is_match * $ranks->index($safe_idx) + ( !$is_match ) * $interp;
    }

    my $out = $out_t->transpose->sever;
    return $is_1d ? $out->slice(':,(0)')->sever : $out;
}

sub prepare_knn_model {
    my ( $class, $w, $u, $t ) = @_;
    require WoLFPSort::UniformScaler;
    my $raw_X = zeroes( $t->{n_feat}, $t->{n_rows} );
    ${ $raw_X->get_dataref } = $t->{packed_buf};
    $raw_X->upd_data;
    my $scaler = WoLFPSort::UniformScaler->new($raw_X);
    return {
        w1        => double( $w->{w1} ),
        w2        => double( $w->{w2} ),
        utility   => double( $u->{matrix} ),
        class_ids => long( $t->{class_ids} ),
        scaler    => $scaler,
        scaled_X  => $scaler->transform($raw_X),
    };
}

sub predict_knn {
    my ( $class, $args ) = @_;
    my ( $state, $query_vec, $ids, $k, $n_cls, $query_id, $ignore_self )
        = @{$args}{qw(state query_vec ids k n_cls query_id ignore_self)};
    my $q_pdl
        = blessed($query_vec) && $query_vec->isa('PDL')
        ? $query_vec
        : double($query_vec);

    my $q_scaled = $state->{scaler}->transform($q_pdl);
    my $diff     = $state->{scaled_X} - $q_scaled;

    # Horton et al. (APBC06), Section 2.3.2 (Classification Algorithm), p. 4,
    # para. 1, Equation 1:
    # d_W(j, k) = sum_i w1_i * |F_ji - F_ki| + sqrt(sum_i w2_i * (F_ji - F_ki)^2)
    my $l1   = sumover( $state->{w1} * abs $diff );
    my $l2   = sqrt sumover( $state->{w2} * ( $diff * $diff ) );
    my $dist = $l1 + $l2;

    my $order = qsorti($dist);
    my $n     = scalar @{$ids};
    $k = $n if $k > $n;

    my @neighbor_indices;
    for my $pos ( 0 .. $n - 1 ) {
        my $idx = $order->at($pos);
        if ( $ignore_self && defined $query_id && $ids->[$idx] eq $query_id ) {
            next;
        }
        push @neighbor_indices, $idx;
        last if @neighbor_indices == $k;
    }

    my $counts = zeroes($n_cls);
    for my $idx (@neighbor_indices) {
        my $c = $state->{class_ids}->at($idx);
        $counts->set( $c, $counts->at($c) + 1.0 );
    }

    # Horton et al. (APBC06), Section 2.3.3 (Extensions for Dual Localization
    # Prediction), p. 4, para. 1 & Table 3: symmetric utility matrix scoring.
    my @scores   = sumover( $state->{utility} * $counts )->list;
    my @nb_dists = map { $dist->at($_) } @neighbor_indices;
    return ( \@scores, \@neighbor_indices, \@nb_dists );
}

sub vector_sum {
    my ( $class, $vals ) = @_;
    return sumover( double($vals) )->sclr;
}

sub weighted_average {
    my ( $class, $weights, $vals ) = @_;
    my $w = double($weights);
    my $v = double($vals);
    return sumover( $w * $v )->sclr / sumover($w)->sclr;
}

sub max_window_sum {
    my ( $class, $vals, $window_len, $n_win ) = @_;
    my $vec = double($vals);
    return sumover($vec)->sclr if $n_win <= 1;
    my $windows = $vec->lags( 0, 1, $window_len );
    return max( sumover( $windows->transpose ) )->sclr;
}

sub encode_indices {
    my ( $class, $idx_ref ) = @_;
    return long($idx_ref);
}

sub psortii_prepare {
    my ( $class, $matrix ) = @_;
    my $mat
        = blessed($matrix) && $matrix->isa('PDL')
        ? $matrix
        : double($matrix);
    croak 'Expected 2D training matrix (F x N)' unless $mat->ndims == 2;
    return $mat;
}

sub psortii_distances {
    my ( $class, $mat, $query ) = @_;
    my $q
        = blessed($query) && $query->isa('PDL')
        ? $query
        : double($query);
    my $diff = $mat - $q;
    return [ sumover( $diff * $diff )->list ];
}

1;

__END__

=head1 NAME

WoLFPSort::Backend::PDL - PDL-vectorized numerical backend for WoLFPSort

=head1 TYPES

Operates on real-valued C<double> ndarrays and C<long> index vectors. Complex
PDL types are not supported.

=cut
