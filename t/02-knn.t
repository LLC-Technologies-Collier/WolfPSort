#!/usr/bin/env perl
# Copyright (C) 2004-2006 Paul B. Horton
# Copyright 2026 Google LLC
# Licensed under the GNU General Public License v3 (GPLv3).
use strict;
use warnings;
use Carp qw(croak);
use Test::More;
use WoLFPSort::KNN;

subtest 'fungi model parity with baseline wolfPredict' => sub {
    my $knn = WoLFPSort::KNN->new(
        weight_file   => 'data/fungi.wolfw',
        training_file => 'data/fungi.wolff',
        utility_file  => 'data/fungi.wolfu',
    );

    is( $knn->k, 27, 'fungi k is 27' );
    ok( $knn->n_samples > 2000, 'fungi training set loaded' );

    open my $fh, '<', 'data/fungi.wolff' or croak "Cannot open data/fungi.wolff: $!";
    my @lines = <$fh>;
    close $fh or croak "Cannot close data/fungi.wolff: $!";

    my @rcc1;
    for my $line (@lines) {
        next unless $line =~ /^RCC1_YEAST\s+/ax;
        $line =~ s/\#.*//ax;
        my @f = split /\s+/ax, $line;
        splice @f, 0, 2;
        @rcc1 = @f;
        last;
    }

    my $res = $knn->predict_vector( \@rcc1, id => 'RCC1_YEAST' );
    is( $res->{scores}[0]{class}, 'nucl', 'RCC1_YEAST top class is nucl' );
    is( WoLFPSort::KNN::format_number( $res->{scores}[0]{score} ),
        '12.5', 'RCC1_YEAST nucl score is 12.5',
    );
    is( $res->{neighbors}[0]{id}, 'RCC1_YEAST', 'nearest neighbor is self' );
    is( $res->{neighbors}[1]{id}, 'TREA_YEAST', '2nd neighbor is TREA_YEAST' );
    is( WoLFPSort::KNN::format_number( $res->{neighbors}[1]{similarity} ),
        '-349.235', 'exact similarity to TREA_YEAST',
    );

    my $res_no_self = $knn->predict_vector( \@rcc1, id => 'RCC1_YEAST', ignore_self => 1 );
    is( $res_no_self->{neighbors}[0]{id}, 'TREA_YEAST', 'ignore_self skips RCC1_YEAST' );
};

subtest 'error paths' => sub {
    my $ok1 = eval { WoLFPSort::KNN->new( weight_file => 'data/fungi.wolfw' ); 1 };
    ok( !$ok1, 'new() fails on missing args' );
    like( $@, qr/Missing required argument/a, 'croaks on missing args' );

    my $ok2 = eval {
        WoLFPSort::KNN->new(
            weight_file   => 'nonexistent.wolfw',
            training_file => 'data/fungi.wolff',
            utility_file  => 'data/fungi.wolfu',
        );
        1;
    };
    ok( !$ok2, 'new() fails on missing weight file' );
    like( $@, qr/Cannot open weight file/a, 'croaks on missing weight file' );
};

done_testing();
