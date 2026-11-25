#!/usr/bin/env perl
# Copyright (C) 2004-2006 Paul B. Horton
# Copyright 2026 Google LLC
# Licensed under the GNU General Public License v3 (GPLv3).
use strict;
use warnings;
use Test::More;
use WoLFPSort::PSORTII::KNN;

subtest 'Euclidean kNN prediction and tie breaking' => sub {
    my $knn = WoLFPSort::PSORTII::KNN->new(
        matrix    => [ [ 0.0, 0.0 ], [ 0.1, 0.1 ], [ 5.0, 5.0 ], [ 5.1, 5.1 ], ],
        names     => [qw( s1 s2 s3 s4 )],
        locations => [qw( nuc nuc mit mit )],
    );

    my $res = $knn->predict( [ 0.05, 0.05 ], 2 );
    is( $res->{prediction},          'nuc', 'predicts nuc for query near (0,0)' );
    is( $res->{classes}[0]{percent}, 100.0, '100% nuc' );

    my $res2 = $knn->predict( [ 4.9, 5.0 ], 2 );
    is( $res2->{prediction}, 'mit', 'predicts mit for query near (5,5)' );
};

done_testing();
