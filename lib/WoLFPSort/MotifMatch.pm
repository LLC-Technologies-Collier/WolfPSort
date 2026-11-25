# Copyright (C) 2004-2006 Paul B. Horton
# Copyright 2026 Google LLC
# Licensed under the GNU General Public License v3 (GPLv3).

package WoLFPSort::MotifMatch;

use strict;
use warnings;

our $VERSION = '0.200';

sub new {
    my ( $class, $description, $offset, $length ) = @_;
    return bless {
        description => $description,
        offset      => $offset,
        length      => $length,
        },
        ref($class) || $class;
}

sub description {
    my ($self) = @_;
    return $self->{description};
}

sub offset {
    my ($self) = @_;
    return $self->{offset};
}

sub length {    ## no critic (Subroutines::ProhibitBuiltinHomonyms)
    my ($self) = @_;
    return $self->{length};
}

1;

__END__

=head1 NAME

WoLFPSort::MotifMatch - Represents a regex motif match on a protein sequence

=cut
