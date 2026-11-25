#!/usr/bin/env perl
# Copyright 2026 Google LLC
# Licensed under the GNU General Public License v3 (GPLv3).
use strict;
use warnings;
use Carp       qw(croak);
use IPC::Open3 qw(open3);
use Symbol     qw(gensym);
use Test::More;

sub _run_summary {
    my ($org) = @_;
    open my $in_fh, '<', 'bin/testQuery.fasta' or croak "Cannot open testQuery.fasta: $!";
    my $fasta = do { local $/ = undef; <$in_fh> };
    close $in_fh or croak "Cannot close testQuery.fasta: $!";

    my $err = gensym;
    my $pid = open3( my $wtr, my $rdr, $err, $^X, 'bin/runWolfPsortSummary', $org );
    print {$wtr} $fasta;
    close $wtr or croak "Cannot close child stdin: $!";
    my $out = do { local $/ = undef; <$rdr> };
    close $rdr or croak "Cannot close child stdout: $!";
    waitpid $pid, 0;
    return ( $?, $out );
}

subtest 'runWolfPsortSummary end-to-end across fungi, animal, plant' => sub {
    my %expected = (
        fungi => {
            k          => 27,
            RCC1_YEAST => 'nucl 12.5',
            RCL1_YEAST => 'mito 15.5',
            RT04_YEAST => 'mito 16',
            EF1A_ASHGO => 'cyto 27',
        },
        animal => {
            k          => 32,
            RCC1_YEAST => 'nucl 30.5',
            RCL1_YEAST => 'mito 22.5',
            RT04_YEAST => 'cyto 15',
            EF1A_ASHGO => 'cyto 27',
        },
        plant => {
            k          => 14,
            RCC1_YEAST => 'chlo 4',
            RCL1_YEAST => 'chlo 7',
            RT04_YEAST => 'mito 7',
            EF1A_ASHGO => 'cyto 14',
        },
    );

    for my $org (qw(fungi animal plant)) {
        my ( $status, $out ) = _run_summary($org);
        is( $status, 0, "$org exited cleanly" );
        my $k_val = $expected{$org}{k};
        like(
            $out,
            qr/^ [#] \s k \s used \s for \s kNN \s is: \s $k_val /amx,
            "$org reports k=$k_val",
        );
        for my $seq (qw(RCC1_YEAST RCL1_YEAST RT04_YEAST EF1A_ASHGO)) {
            my $exp = $expected{$org}{$seq};
            like( $out, qr/^ $seq \s+ \Q$exp\E /amx, "$org $seq matches baseline top prediction", );
        }
    }
};

done_testing();
