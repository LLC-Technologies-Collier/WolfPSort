#!/usr/bin/env perl
# Copyright (C) 2004-2006 Paul B. Horton
# Copyright 2026 Google LLC
# Licensed under the GNU General Public License v3 (GPLv3).
use strict;
use warnings;
use Carp qw(croak);
use Test::More;
use WoLFPSort::AminoAcidIndexer;
use WoLFPSort::MotifMatch;
use WoLFPSort::RegexMotif;
use WoLFPSort::RegexMotifs;

subtest 'AminoAcidIndexer scalar and PDL encoding' => sub {
    is( WoLFPSort::AminoAcidIndexer::size(),        23, '23 amino acid codes' );
    is( WoLFPSort::AminoAcidIndexer::to_index('A'), 0,  'A -> 0' );
    is( WoLFPSort::AminoAcidIndexer::to_index('X'), 22, 'X -> 22' );

    my $pdl = WoLFPSort::AminoAcidIndexer::encode_sequence('ARNX');
    is_deeply( [ $pdl->list ], [ 0, 1, 2, 22 ], 'encode_sequence returns PDL long ndarray' );

    my $ok = eval { WoLFPSort::AminoAcidIndexer::to_index('?'); 1 };
    ok( !$ok, 'to_index fails on unknown residue' );
    like( $@, qr/Unrecognized \s one-letter \s amino \s acid \s code/ax,
        'rejects unknown residue' );
};

subtest 'RegexMotif and RegexMotifs on real motif library' => sub {
    my $m = WoLFPSort::RegexMotif->new( 'test motif', 'K[RK]{2}' );
    is( $m->num_matches('AAAKRRBBBKKKCCC'), 2, 'finds 2 non-overlapping matches' );

    my @matches = $m->matches('AAAKRRBBBKKKCCC');
    is( scalar @matches,     2, 'returns 2 MotifMatch objects' );
    is( $matches[0]->offset, 3, 'first match offset is 3' );
    is( $matches[0]->length, 3, 'first match length is 3' );

    my $lib_path = 'bin/psortModifiedForWolfFiles/motifLibraries/ribosomal.regexMotifs';
    open my $fh, '<', $lib_path or croak "Cannot open $lib_path: $!";
    my $lib = WoLFPSort::RegexMotifs->new($fh);
    close $fh or croak "Cannot close $lib_path: $!";

    ok( length( $lib->name ) > 0,            'loaded ribosomal motif library name' );
    ok( scalar( @{ $lib->motifs_ref } ) > 0, 'loaded ribosomal motifs' );
};

done_testing();
