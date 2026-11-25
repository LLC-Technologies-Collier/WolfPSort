#!/usr/bin/env perl
# Copyright (C) 2004-2006 Paul B. Horton
# Copyright 2026 Google LLC
# Licensed under the GNU General Public License v3 (GPLv3).
use strict;
use warnings;
use Test::More;

unless ( $ENV{AUTHOR_TESTING} || $ENV{RELEASE_TESTING} ) {
    plan( skip_all => 'Author/Release tests not required for installation' );
}

my $ok = eval { require Test::Kwalitee; Test::Kwalitee->import(); 1 };
plan( skip_all => 'Test::Kwalitee required for testing distribution quality' ) unless $ok;
