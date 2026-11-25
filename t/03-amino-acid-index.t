#!/usr/bin/env perl
# Copyright (C) 2004-2006 Paul B. Horton
# Copyright 2026 Google LLC
# Licensed under the GNU General Public License v3 (GPLv3).
use strict;
use warnings;
use Test::More;
use WoLFPSort::AminoAcidIndex;

subtest 'Kyte-Doolittle and Klein net charge with PDL' => sub {
    my $kd = WoLFPSort::AminoAcidIndex->new('KYTJ820101');
    is( $kd->name,       'Kyte Doolittle', 'index name' );
    is( $kd->value('I'), 4.5,              'I hydropathy is 4.5' );

    my @seq = split //a, 'MKRISTTITTTITITTGNGAGIIIIIIIIII';
    my $sum = $kd->seqAverage( \@seq, 0, 5 );
    ok( abs( $sum - ( 1.9 - 3.9 - 4.5 + 4.5 - 0.8 ) ) < 1e-9, 'seqAverage matches sum' );

    my $max_win = $kd->maxWindowAverage( \@seq, 5, 0, 30 );
    ok( abs( $max_win - ( 5 * 4.5 ) ) < 1e-9, 'maxWindowAverage finds I-rich window' );

    my $klein = WoLFPSort::AminoAcidIndex->new('KLEP840101');
    my $chg   = $klein->seqAverage( \@seq, 0, 5 );
    is( $chg, 2, 'Klein net charge for MKRIS is 2' );
};

subtest 'error paths' => sub {
    my $ok1 = eval { WoLFPSort::AminoAcidIndex->new('NO_SUCH_INDEX'); 1 };
    ok( !$ok1, 'new() fails on unknown index' );
    like( $@, qr/No amino acid index/a, 'rejects unknown index' );

    my $kd  = WoLFPSort::AminoAcidIndex->new('KYTJ820101');
    my $ok2 = eval { $kd->value('?'); 1 };
    ok( !$ok2, 'value() fails on invalid residue' );
    like( $@, qr/not a valid amino acid residue/a, 'rejects invalid residue' );
};

done_testing();
