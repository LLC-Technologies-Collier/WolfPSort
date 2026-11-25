#!/usr/bin/env perl
# Copyright (C) 2004-2006 Paul B. Horton
# Copyright 2026 Google LLC
# Licensed under the GNU General Public License v3 (GPLv3).
use strict;
use warnings;
use Test::More;
use PDL::Lite;
use PDL::Core qw(pdl double);
use WoLFPSort::UniformScaler;

subtest 'basic quantile scaling' => sub {
    my $X      = double( [ [ 10, 100 ], [ 20, 200 ], [ 30, 300 ], [ 40, 400 ] ] );
    my $scaler = WoLFPSort::UniformScaler->new($X);

    is( $scaler->n_features, 2, 'n_features is 2' );
    is( $scaler->n_samples,  4, 'n_samples is 4' );

    my $U = $scaler->transform($X);
    is_deeply( [ $U->dims ], [ 2, 4 ], 'scaled matrix preserves (F, N) shape' );

    my @f0 = $U->slice('(0),:')->list;
    is( $f0[0], 0,  'min value maps to 0' );
    is( $f0[1], 25, '2nd value maps to 25' );
    is( $f0[2], 50, '3rd value maps to 50' );
    is( $f0[3], 75, 'max value maps to 75' );

    my $q = $scaler->transform( double( [ 20, 300 ] ) );
    is_deeply( [ $q->dims ], [2], '1D query returns 1D ndarray' );
    is( $q->at(0), 25, '1D query feature 0' );
    is( $q->at(1), 50, '1D query feature 1' );
};

subtest 'error paths' => sub {
    my $ok1 = eval { WoLFPSort::UniformScaler->new('not_a_pdl'); 1 };
    ok( !$ok1, 'new() fails on non-PDL input' );
    like( $@, qr/Expected a 2D PDL ndarray/a, 'rejects non-PDL input' );

    my $ok2 = eval { WoLFPSort::UniformScaler->new( double( [ 1, 2, 3 ] ) ); 1 };
    ok( !$ok2, 'new() fails on 1D ndarray' );
    like( $@, qr/Expected a 2D PDL ndarray/a, 'rejects 1D ndarray' );

    my $scaler = WoLFPSort::UniformScaler->new( double( [ [ 1, 2 ], [ 3, 4 ] ] ) );
    my $ok3    = eval { $scaler->transform( double( [ 1, 2, 3 ] ) ); 1 };
    ok( !$ok3, 'transform() fails on mismatched features' );
    like( $@, qr/Feature count mismatch/a, 'rejects mismatched feature count' );
};

done_testing();
