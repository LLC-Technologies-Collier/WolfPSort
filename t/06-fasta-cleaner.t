#!/usr/bin/env perl
# Copyright (C) 2004-2006 Paul B. Horton
# Copyright 2026 Google LLC
# Licensed under the GNU General Public License v3 (GPLv3).
use strict;
use warnings;
use Carp qw(croak);
use Test::More;
use WoLFPSort::FastaCleaner;

subtest 'cleans valid and invalid FASTA records' => sub {
    my $input = join "\n",
        '>seq1 valid protein',
        'mkristtitttitittgngagiiiiiiiiiiaaaaa',
        '>seq2 too short',
        'MKRIST',
        '>seq3 has pyrrolysine O',
        'MKRISTTITTTITITTGNGAGIIIIIIIIIIOAAAA',
        '';

    my $output  = '';
    my $cleaner = WoLFPSort::FastaCleaner->new();
    my $ok;
    {
        open my $in_fh,  '<', \$input  or croak "Cannot open input scalar: $!";
        open my $out_fh, '>', \$output or croak "Cannot open output scalar: $!";
        $ok = $cleaner->clean_fh( $in_fh, $out_fh );
        close $in_fh  or croak "Cannot close input scalar: $!";
        close $out_fh or croak "Cannot close output scalar: $!";
    }

    ok( $ok, 'clean_fh returns true' );
    like( $output, qr/^ >seq1 \s valid \s protein \n MKRIST /amx, 'uppercases valid sequence' );
    like( $output, qr/Sequence \s seq2 \s skipped \s because \s length \s 6/ax, 'warns on short' );
    like( $output, qr/Sequence \s seq3 \s skipped \s because \s it \s contains/ax, 'warns on O' );
};

done_testing();
