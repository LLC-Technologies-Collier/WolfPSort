#!/usr/bin/env perl
# Copyright (C) 2004-2006 Paul B. Horton
# Copyright 2026 Google LLC
# Licensed under the GNU General Public License v3 (GPLv3).
use strict;
use warnings;
use WoLFPSort::CLI::Predict;

exit WoLFPSort::CLI::Predict::run( \@ARGV, \*STDIN, \*STDOUT );
