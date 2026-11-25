#!/usr/bin/env perl
# Copyright (C) 2004-2006 Paul B. Horton
# Copyright 2026 Google LLC
# Licensed under the GNU General Public License v3 (GPLv3).
use strict;
use warnings;
use Carp qw(croak);
use Test::More;

my $has_xs = eval { require WoLFPSort::XS; 1 };
plan( skip_all => 'WoLFPSort::XS not built' ) unless $has_xs;

subtest 'WoLFPSort::XS in-process C99 parity on fungi' => sub {
    my $m = WoLFPSort::XS->new( 'data/fungi.wolfw', 'data/fungi.wolff', 'data/fungi.wolfu', );
    is( $m->k, 27, 'XS fungi k is 27' );
    ok( $m->n_samples > 2000, 'XS fungi training set loaded' );

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

    my $res = $m->predict_vector( \@rcc1, 'RCC1_YEAST', 0 );
    is( $res->{scores}[0]{class},                    'nucl', 'XS RCC1_YEAST top class is nucl' );
    is( sprintf( '%.1f', $res->{scores}[0]{score} ), '12.5', 'XS nucl score is 12.5' );
    is( $res->{neighbors}[0]{id}, 'RCC1_YEAST',              'XS nearest neighbor is self' );
    is( $res->{neighbors}[1]{id}, 'TREA_YEAST',              'XS 2nd neighbor is TREA_YEAST' );
    is( sprintf( '%.3f', $res->{neighbors}[1]{similarity} ), '-349.235', 'XS exact similarity' );

    my $res_ns = $m->predict_vector( \@rcc1, 'RCC1_YEAST', 1 );
    is( $res_ns->{neighbors}[0]{id}, 'TREA_YEAST', 'XS ignore_self skips RCC1_YEAST' );

    my $ok = eval { $m->predict_vector( [ 1, 2, 3 ], 'bad', 0 ); 1 };
    ok( !$ok, 'XS croaks on feature count mismatch' );
    like( $@, qr/Feature \s count \s mismatch/ax, 'XS mismatch error message' );
};

done_testing();
