# Copyright (C) 2004-2006 Paul B. Horton
# Copyright 2026 Google LLC
# Licensed under the GNU General Public License v3 (GPLv3).

package WoLFPSort::RegexMotifs;

use strict;
use warnings;
use Carp qw(croak);
use WoLFPSort::RegexMotif;

our $VERSION = '0.200';

sub new {
    my ( $class, $istream )    = @_;
    my ( $name,  $motifs_ref ) = read_motifs_from_stream($istream);
    return bless {
        name       => $name,
        motifs_ref => $motifs_ref,
        },
        ref($class) || $class;
}

sub name {
    my ($self) = @_;
    return $self->{name};
}

sub motifs_ref {
    my ($self) = @_;
    return $self->{motifs_ref};
}

sub motifsRef {
    my ($self) = @_;
    return $self->{motifs_ref};
}

sub read_motifs_from_stream {
    my ($istream) = @_;
    my @motifs;
    my ( $description, $pattern, $cur_line );
    my $num_lines_read = 0;
    my $name           = '';

    while ( my $line = <$istream> ) {
        chomp $line;
        $cur_line = $line;
        next if $cur_line     =~ /^\#/ax;
        next unless $cur_line =~ /\S+/ax;

        unless ($name) {
            $name = $line;
            next;
        }

        if ( $num_lines_read == 0 ) {
            $description = $line;
        }
        else {
            $pattern = $line;
            push @motifs, WoLFPSort::RegexMotif->new( $description, $pattern );
        }
        $num_lines_read = ( ++$num_lines_read % 2 );
    }

    croak "Could not parse motif line: $cur_line" if $num_lines_read;
    return ( $name, \@motifs );
}

sub num_matches {
    my ( $self, $seq ) = @_;
    my $total = 0;
    for my $motif ( @{ $self->motifs_ref() } ) {
        $total += $motif->num_matches($seq);
    }
    return $total;
}

sub numMatches {
    my ( $self, $seq ) = @_;
    return $self->num_matches($seq);
}

sub report_matches {
    my ( $self, $seq, $endl ) = @_;
    $endl ||= "\n";

    my $total = 0;
    my $mots  = $self->motifs_ref();
    my $name  = $self->name();
    print 'checking ', scalar(@$mots), " $name: ";

    for my $motif (@$mots) {
        my @matches = $motif->matches($seq);
        next unless @matches;
        $total += scalar @matches;

        my $desc = $motif->description();
        print "$endl      $desc:  *** found ***$endl";
        for my $match (@matches) {
            my $sub = substr $seq, $match->offset(), $match->length();
            printf "         %s at %d$endl", $sub, $match->offset() + 1;
        }
    }

    print " none$endl" unless $total;
    print $endl;
    return $total;
}

sub reportMatches {
    my ( $self, $seq, $endl ) = @_;
    return $self->report_matches( $seq, $endl );
}

1;

__END__

=head1 NAME

WoLFPSort::RegexMotifs - Collection of regular-expression protein sequence motifs

=cut
