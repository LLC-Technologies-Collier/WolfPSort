#!/usr/bin/env perl
# Copyright (C) 2004-2006 Paul B. Horton
# Copyright 2026 Google LLC
# Licensed under the GNU General Public License v3 (GPLv3).
use strict;
use warnings;

use WoLFPSort::FastaCleaner;

my $usage = "$0 [--html]";
my $html  = 0;
my $arg   = shift;

if ( defined $arg ) {
    if ( $arg eq '--html' ) {
        $html = 1;
    }
    else {
        print "Error, unparsed arg: $arg\n\n$usage\n";
        exit 255;
    }
}

my $cleaner = WoLFPSort::FastaCleaner->new( html => $html );
$cleaner->clean_fh( \*STDIN, \*STDOUT );
exit 0;
