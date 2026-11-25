# Copyright (C) 2004-2006 Paul B. Horton
# Copyright 2026 Google LLC
# Licensed under the GNU General Public License v3 (GPLv3).

package WoLFPSort::AminoAcidIndex;

use strict;
use warnings;
use Carp qw(croak);
use WoLFPSort::Backend;

our $VERSION = '0.200';

my %relativeAbundance = (
    name => 'relative abundance',
    L    => 0.0994,
    A    => 0.0884,
    G    => 0.0703,
    V    => 0.0677,
    S    => 0.0672,
    E    => 0.0624,
    I    => 0.0595,
    R    => 0.0570,
    T    => 0.0543,
    D    => 0.0539,
    K    => 0.0527,
    P    => 0.0471,
    N    => 0.0417,
    F    => 0.0400,
    Q    => 0.0382,
    Y    => 0.0300,
    M    => 0.0237,
    H    => 0.0220,
    C    => 0.0124,
    W    => 0.0121,
);

my %EngelmanGES1986 = (
    name => 'Engelman GES 1986',
    A    => -6.7,
    R    =>  51.5,
    N    =>  20.1,
    D    =>  38.5,
    C    => -8.4,
    Q    =>  17.2,
    E    =>  34.3,
    G    => -4.2,
    H    =>  12.6,
    I    => -13.0,
    L    => -11.7,
    K    =>  36.8,
    M    => -14.2,
    F    => -15.5,
    P    =>  0.8,
    S    => -2.5,
    T    => -5.0,
    W    => -7.9,
    Y    =>  2.9,
    V    => -10.9,
    B    =>  30.5,
    Z    =>  27.8,
    X    =>  5.8,
);

my %KYTJ820101 = (
    name => 'Kyte Doolittle',
    id   => 'KYTJ820101',
    A    =>  1.8,
    R    => -4.5,
    N    => -3.5,
    D    => -3.5,
    C    =>  2.5,
    Q    => -3.5,
    E    => -3.5,
    G    => -0.4,
    H    => -3.2,
    I    =>  4.5,
    L    =>  3.8,
    K    => -3.9,
    M    =>  1.9,
    F    =>  2.8,
    P    => -1.6,
    S    => -0.8,
    T    => -0.7,
    W    => -0.9,
    Y    => -1.3,
    V    =>  4.2,
    B    => -3.5,
    Z    => -3.5,
    X    => -0.19,
);

my %KLEP840101 = (
    name => 'Klein net Charge',
    id   => 'KLEP840101',
    A    =>  0,
    R    =>  1,
    N    =>  0,
    D    => -1,
    C    =>  0,
    Q    =>  0,
    E    => -1,
    G    =>  0,
    H    =>  0,
    I    =>  0,
    L    =>  0,
    K    =>  1,
    M    =>  0,
    F    =>  0,
    P    =>  0,
    S    =>  0,
    T    =>  0,
    W    =>  0,
    Y    =>  0,
    V    =>  0,
    B    => -0.56,
    Z    => -0.62,
    X    => -0.01,
);

my %nameToReference = (
    engelmanges1986 => \%EngelmanGES1986,
    kytedoolittle   => \%KYTJ820101,
    kleinnetcharge  => \%KLEP840101,
);

my %idToReference = (
    engelmanges1986 => \%EngelmanGES1986,
    KLEP840101      => \%KLEP840101,
    KYTJ820101      => \%KYTJ820101,
);

my @standard20Residues = qw( A C D E F G H I K L M N P Q R S T V W Y );

sub new {
    my ( $class, $indexIdOrName ) = @_;
    my $ref = _idToReference($indexIdOrName) || _nameToReference($indexIdOrName);
    croak "No amino acid index with id or name: '$indexIdOrName' known" unless $ref;
    return bless { %{$ref} }, $class;
}

sub name {
    my ($self) = @_;
    return $self->{name};
}

sub id {
    my ($self) = @_;
    return $self->{id};
}

sub value {
    my ( $self, $residue, @extra ) = @_;
    croak 'expected 2 arguments' if @extra || !defined $residue;
    croak "'$residue' not a valid amino acid residue"
        unless length($residue) == 1 && exists $self->{$residue};
    return $self->{$residue};
}

sub indexNames {
    return keys %nameToReference;
}

sub abundanceWeightedAverageValue {
    my ( $self, @extra ) = @_;
    croak 'expected one argument' if @extra;
    my @abund = map { $relativeAbundance{$_} } @standard20Residues;
    my @vals  = map { $self->value($_) } @standard20Residues;
    return WoLFPSort::Backend->current()->weighted_average( \@abund, \@vals );
}

sub seqAverage {
    my ( $self, @args ) = @_;
    my $seqArrayRef = $args[0];
    my ( $startPos, $endPos ) = _startEndPos(@args);
    return 0 if $endPos <= $startPos;
    my @vals = map { $self->value( $seqArrayRef->[$_] ) } ( $startPos .. $endPos - 1 );
    return WoLFPSort::Backend->current()->vector_sum( \@vals );
}

sub maxWindowAverage {
    my ( $self, $seqArrayRef, $windowLength, @rest ) = @_;
    my ( $startPos, $endPos ) = _startEndPos( $seqArrayRef, @rest );

    my $numWindowPositions = $endPos - $startPos - $windowLength + 1;
    croak "startPos: $startPos, endPos: $endPos, windowLength: $windowLength"
        if $numWindowPositions < 0;

    my $last_start = $numWindowPositions > $startPos ? $numWindowPositions : $startPos;
    my $slice_end  = $last_start + $windowLength;
    my $n_win      = $last_start - $startPos + 1;
    my @vals       = map { $self->value( $seqArrayRef->[$_] ) } ( $startPos .. $slice_end - 1 );

    return WoLFPSort::Backend->current()->max_window_sum( \@vals, $windowLength, $n_win );
}

sub _startEndPos {
    my (@args) = @_;
    croak 'expected arguments' unless @args;
    my $seqArrayRef = shift @args;
    push @args, 0                      if @args == 0;
    push @args, scalar @{$seqArrayRef} if @args == 1;
    croak 'expected 1 to 3 arguments' unless @args == 2;
    my ( $startPos, $length ) = @args;
    my $seqLen = scalar @{$seqArrayRef};

    croak 'invalid startPos value' if $startPos < 0 || $startPos > $seqLen;
    my $endPos = $length < 0 ? $seqLen - $length : $startPos + $length;
    croak "endPos: $endPos > seqLen: $seqLen" if $endPos > $seqLen;
    return ( $startPos, $endPos );
}

sub _nameToReference {
    my ($indexName) = @_;
    return 0 unless defined $indexName;
    $indexName =~ tr/ //d;
    $indexName =~ tr/A-Z/a-z/;
    return exists $nameToReference{$indexName} ? $nameToReference{$indexName} : 0;
}

sub _idToReference {
    my ($id) = @_;
    return 0 unless defined $id;
    return exists $idToReference{$id} ? $idToReference{$id} : 0;
}

1;

__END__

=head1 NAME

WoLFPSort::AminoAcidIndex - Amino acid index and sliding-window sums

=head1 SYNOPSIS

    use WoLFPSort::AminoAcidIndex;

    my $kd = WoLFPSort::AminoAcidIndex->new('KYTJ820101');
    my $max_hydro = $kd->maxWindowAverage(\@residues, 12, 0, 29);

=cut
