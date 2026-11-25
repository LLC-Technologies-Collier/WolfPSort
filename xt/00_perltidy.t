#!/usr/bin/env perl
# Copyright (C) 2004-2006 Paul B. Horton
# Copyright 2026 Google LLC
# Licensed under the GNU General Public License v3 (GPLv3).
use strict;
use warnings;
use Test::More;

unless ( $ENV{TEST_PERLTIDY} ) {
    plan( skip_all => 'Set TEST_PERLTIDY=1 to enable perltidy tests' );
}

my $ok = eval { require Test::PerlTidy; Test::PerlTidy->import(); 1 };
plan( skip_all => 'Test::PerlTidy required for checking code formatting' ) unless $ok;

run_tests( exclude => [ qr{^bin/(?:fastafmt|psortModifiedForWolfFiles)/}ax, qr{^tmp/}ax ] );
