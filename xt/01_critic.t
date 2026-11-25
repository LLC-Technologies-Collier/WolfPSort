#!/usr/bin/env perl
# Copyright 2026 Google LLC
# Licensed under the GNU General Public License v3 (GPLv3).

use strict;
use warnings;
use Test::More;
use Test::Perl::Critic -profile => '.perlcriticrc', -exclude => ['ProhibitUselessNoCritic'];
use File::Find;

my @files;
find(
    sub {
        push @files, $File::Find::name if /\.pm$/ && -f $_;
    },
    'lib'
);
push @files, 'bin/runWolfPsortSummary', 'bin/wolfPredict.pl', 'bin/checkFastaInput.pl',
    glob 't/*.t';

plan tests => scalar @files;

foreach my $file (@files) {
    critic_ok( $file, "$file is critic clean" );
}
