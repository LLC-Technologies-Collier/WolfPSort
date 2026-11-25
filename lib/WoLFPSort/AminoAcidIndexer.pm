# Copyright (C) 2004-2006 Paul B. Horton
# Copyright 2026 Google LLC
# Licensed under the GNU General Public License v3 (GPLv3).

package WoLFPSort::AminoAcidIndexer;

use strict;
use warnings;
use Carp qw(croak);
use WoLFPSort::Backend;

our $VERSION = '0.200';

my %CODE_TO_INDEX = (
    'A' => 0,
    'R' => 1,
    'N' => 2,
    'D' => 3,
    'C' => 4,
    'Q' => 5,
    'E' => 6,
    'G' => 7,
    'H' => 8,
    'I' => 9,
    'L' => 10,
    'K' => 11,
    'M' => 12,
    'F' => 13,
    'P' => 14,
    'S' => 15,
    'T' => 16,
    'W' => 17,
    'Y' => 18,
    'V' => 19,
    'B' => 20,
    'Z' => 21,
    'X' => 22,
);

sub to_index {
    my ($code) = @_;
    croak 'to_index expects one amino acid code argument'
        unless defined $code && length $code;
    return $CODE_TO_INDEX{$code} if exists $CODE_TO_INDEX{$code};
    croak "Unrecognized one-letter amino acid code: '$code'";
}

sub toIndex {
    my ($code) = @_;
    return to_index($code);
}

sub size {
    return scalar keys %CODE_TO_INDEX;
}

sub encode_sequence {
    my ($seq) = @_;
    croak 'Expected a non-empty sequence string' unless defined $seq && length $seq;
    my @idx = map { to_index($_) } split //a, $seq;
    return WoLFPSort::Backend->current()->encode_indices( \@idx );
}

1;

__END__

=head1 NAME

WoLFPSort::AminoAcidIndexer - Map one-letter amino acid codes to integer indices

=head1 SYNOPSIS

    use WoLFPSort::AminoAcidIndexer;

    my $idx = WoLFPSort::AminoAcidIndexer::to_index('A'); # 0
    my $enc = WoLFPSort::AminoAcidIndexer::encode_sequence('MKRIS');

=cut
