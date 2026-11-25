# Copyright (C) 2004-2006 Paul B. Horton
# Copyright 2026 Google LLC
# Licensed under the GNU General Public License v3 (GPLv3).

package WoLFPSort::RegexMotif;

use strict;
use warnings;
use WoLFPSort::MotifMatch;

our $VERSION = '0.200';

use overload
    q{""} => sub {
    my ($self) = @_;
    my $d      = $self->description();
    my $p      = $self->pattern();
    return "description: $d\npattern: $p\n";
    },
    'bool'   => sub { return 1; },
    fallback => 1;

sub new {
    my ( $class, $description, $pattern ) = @_;
    return bless {
        description => $description,
        pattern     => $pattern,
        },
        ref($class) || $class;
}

sub description {
    my ($self) = @_;
    return $self->{description};
}

sub pattern {
    my ($self) = @_;
    return $self->{pattern};
}

sub matches {
    my ( $self, $seq ) = @_;
    my $pat        = $self->pattern();
    my @match_strs = $seq =~ /$pat/agx;
    my @matches;
    my $search_pos = 0;

    for my $m (@match_strs) {
        my $pos = index $seq, $m, $search_pos;
        $search_pos = $pos + 1;
        push @matches, WoLFPSort::MotifMatch->new( $self->description(), $pos, length $m, );
    }
    return @matches;
}

sub num_matches {
    my ( $self, $seq ) = @_;
    my $pat     = $self->pattern();
    my @matches = $seq =~ /$pat/agx;
    return scalar @matches;
}

sub numMatches {
    my ( $self, $seq ) = @_;
    return $self->num_matches($seq);
}

1;

__END__

=head1 NAME

WoLFPSort::RegexMotif - Regular-expression protein sequence motif

=cut
