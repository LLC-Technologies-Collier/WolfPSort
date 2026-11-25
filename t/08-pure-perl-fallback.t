#!/usr/bin/env perl
# Copyright (C) 2004-2006 Paul B. Horton
# Copyright 2026 Google LLC
# Licensed under the GNU General Public License v3 (GPLv3).
use strict;
use warnings;
use Carp qw(croak);
use Test::More;

use WoLFPSort;
use WoLFPSort::UniformScaler;
use WoLFPSort::KNN;
use WoLFPSort::AminoAcidIndex;
use WoLFPSort::AminoAcidIndexer;
use WoLFPSort::PSORTII::KNN;

local $ENV{WOLFPSORT_NO_PDL} = 1;
local $ENV{WOLFPSORT_NO_XS}  = 1;

ok( !WoLFPSort::has_pdl(), 'WOLFPSORT_NO_PDL=1 disables PDL backend' );

subtest 'pure-Perl UniformScaler' => sub {
    my $rows   = [ [ 10, 100 ], [ 20, 200 ], [ 30, 300 ], [ 40, 400 ] ];
    my $scaler = WoLFPSort::UniformScaler->new($rows);
    my $scaled = $scaler->transform($rows);
    is_deeply(
        [ map { $_->[0] } @{$scaled} ],
        [ 0, 25, 50, 75 ],
        'pure-Perl mid-rank percentiles match PDL',
    );
    my $q = $scaler->transform( [ 20, 300 ] );
    is_deeply( $q, [ 25, 50 ], '1D pure-Perl query transform' );
};

subtest 'pure-Perl KNN parity with fungi baseline' => sub {
    my $knn = WoLFPSort::KNN->new(
        weight_file   => 'data/fungi.wolfw',
        training_file => 'data/fungi.wolff',
        utility_file  => 'data/fungi.wolfu',
    );

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
    is( $res->{scores}[0]{class}, 'nucl', 'pure-Perl top class is nucl' );
    is( WoLFPSort::KNN::format_number( $res->{scores}[0]{score} ),
        '12.5', 'pure-Perl nucl score is 12.5',
    );
    is( $res->{neighbors}[1]{id}, 'TREA_YEAST', 'pure-Perl 2nd neighbor is TREA_YEAST' );
    is( WoLFPSort::KNN::format_number( $res->{neighbors}[1]{similarity} ),
        '-349.235', 'pure-Perl exact similarity to TREA_YEAST',
    );
};

subtest 'pure-Perl AminoAcidIndex, Indexer, and PSORTII::KNN' => sub {
    my $kd  = WoLFPSort::AminoAcidIndex->new('KYTJ820101');
    my @seq = split //a, 'MKRISTTITTTITITTGNGAGIIIIIIIIII';
    ok( abs( $kd->maxWindowAverage( \@seq, 5, 0, 30 ) - 22.5 ) < 1e-9,
        'pure-Perl maxWindowAverage' );

    my $enc = WoLFPSort::AminoAcidIndexer::encode_sequence('ARNX');
    is_deeply( $enc, [ 0, 1, 2, 22 ], 'pure-Perl encode_sequence returns arrayref' );

    my $p2 = WoLFPSort::PSORTII::KNN->new(
        matrix    => [ [ 0.0, 0.0 ], [ 0.1, 0.1 ], [ 5.0, 5.0 ], [ 5.1, 5.1 ] ],
        names     => [qw( s1 s2 s3 s4 )],
        locations => [qw( nuc nuc mit mit )],
    );
    my $res = $p2->predict( [ 0.05, 0.05 ], 2 );
    is( $res->{prediction}, 'nuc', 'pure-Perl PSORTII::KNN predicts nuc' );
};

done_testing();
