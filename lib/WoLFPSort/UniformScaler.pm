# Copyright (C) 2004-2006 Paul B. Horton
# Copyright 2026 Google LLC
# Licensed under the GNU General Public License v3 (GPLv3).

package WoLFPSort::UniformScaler;

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
    require PDL::Slices;
    require PDL::Primitive;
    $_PDL_LOADED = 1;
    return;
}

sub _build_pdl_tables {
    my ( $matrix, $n_features, $n_samples ) = @_;
    _ensure_pdl();
    my $sorted = $matrix->transpose->qsort->sever;
    my @tables;
    for my $f ( 0 .. $n_features - 1 ) {
        my $s_col = $sorted->slice(":,($f)");
        my ( $counts, $vals ) = PDL::Slices::rle($s_col);
        my $nz = PDL::Primitive::which( $counts > 0 );
        my $c  = $counts->index($nz)->double;
        my $v  = $vals->index($nz)->sever;
        my $hi = PDL::Ufunc::cumusumover($c);
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
    return \@tables;
}

sub _build_pp_tables {
    my ( $rows, $n_features, $n_samples ) = @_;
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
    return \@tables;
}

sub new {
    my ( $class, $matrix ) = @_;

    if ( blessed($matrix) && $matrix->isa('PDL') ) {
        croak 'Expected a 2D PDL ndarray (F x N)' unless $matrix->ndims == 2;
        my ( $n_features, $n_samples ) = $matrix->dims;
        croak 'Training matrix must have at least 1 feature and 1 sample'
            if $n_features < 1 || $n_samples < 1;

        return bless {
            backend    => 'pdl',
            tables     => _build_pdl_tables( $matrix, $n_features, $n_samples ),
            n_features => $n_features,
            n_samples  => $n_samples,
        }, $class;
    }

    if ( ref($matrix) eq 'ARRAY' && @{$matrix} && ref( $matrix->[0] ) eq 'ARRAY' ) {
        my $n_samples  = scalar @{$matrix};
        my $n_features = scalar @{ $matrix->[0] };
        croak 'Training matrix must have at least 1 feature and 1 sample'
            if $n_features < 1 || $n_samples < 1;

        return bless {
            backend    => 'pp',
            tables     => _build_pp_tables( $matrix, $n_features, $n_samples ),
            n_features => $n_features,
            n_samples  => $n_samples,
        }, $class;
    }

    croak 'Expected a 2D PDL ndarray (F x N) or arrayref of feature rows';
}

sub n_features {
    my ($self) = @_;
    return $self->{n_features};
}

sub n_samples {
    my ($self) = @_;
    return $self->{n_samples};
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
# sentences 3-4: normalize each feature to its percentile value in the training
# data; for unseen test values, use linear interpolation (and extrapolation
# outside [min, max]) to obtain a pseudo-percentile value.
sub _interp_pp {
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

sub _transform_pp_row {
    my ( $self, $row ) = @_;
    my $n_feat = $self->{n_features};
    croak 'Feature count mismatch: got ' . scalar( @{$row} ) . ", expected $n_feat"
        if @{$row} != $n_feat;

    my $tables = $self->{tables};
    my @out;
    for my $f ( 0 .. $n_feat - 1 ) {
        my $tbl = $tables->[$f];
        my $x   = $row->[$f] + 0;
        my $n_u = $tbl->{n_u};
        my $idx = _bsearch_left( $x, $tbl->{vals}, $n_u );
        if ( $idx < $n_u && $tbl->{vals}[$idx] == $x ) {
            push @out, $tbl->{ranks}[$idx];
        }
        else {
            push @out, _interp_pp( $x, $tbl->{vals}, $tbl->{ranks}, $n_u, $idx );
        }
    }
    return \@out;
}

sub transform {
    my ( $self, $X ) = @_;

    if ( $self->{backend} eq 'pp' ) {
        croak 'Expected an arrayref of features or arrayref of rows'
            unless ref($X) eq 'ARRAY';
        if ( @{$X} && ref( $X->[0] ) eq 'ARRAY' ) {
            return [ map { $self->_transform_pp_row($_) } @{$X} ];
        }
        return $self->_transform_pp_row($X);
    }

    croak 'Expected a 1D or 2D PDL ndarray'
        unless blessed($X) && $X->isa('PDL');

    _ensure_pdl();
    my $is_1d = ( $X->ndims == 1 );
    my $mat   = $is_1d ? $X->dummy(1) : $X;
    croak 'Expected 1D (F) or 2D (F x M) PDL ndarray'
        unless $mat->ndims == 2;

    my ( $n_feat, $n_obs ) = $mat->dims;
    my $expected_f = $self->{n_features};
    croak "Feature count mismatch: got $n_feat, expected $expected_f"
        if $n_feat != $expected_f;

    my $N      = $self->{n_samples};
    my $tables = $self->{tables};
    my $mat_t  = $mat->transpose->sever;
    my $out_t  = PDL::Core::zeroes( $n_obs, $n_feat );

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
        my $idx       = PDL::Primitive::vsearch( $x_col, $v, { mode => 'insert_leftmost' } );
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

1;

__END__

=head1 NAME

WoLFPSort::UniformScaler - Empirical-CDF feature scaler for WoLF PSORT

=head1 SYNOPSIS

    use WoLFPSort::UniformScaler;

    my $scaler = WoLFPSort::UniformScaler->new($train_matrix);
    my $scaled_query = $scaler->transform($query_vector);

=head1 DESCRIPTION

Maps each feature column to uniform C<[0, 100]> mid-rank percentiles, matching
C<cbrc::UniformScaler> in the original C++ C<wolfPredict>. Uses L<PDL> when
given a 2D PDL ndarray and falls back to a pure-Perl RLE + binary search
implementation when given an arrayref of rows.

=head1 TYPES

When using L<PDL>, operates on real-valued C<double> ndarrays. Complex types
are not supported.

=cut
