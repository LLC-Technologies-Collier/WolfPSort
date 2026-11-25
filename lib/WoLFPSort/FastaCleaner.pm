# Copyright (C) 2004-2006 Paul B. Horton
# Copyright 2026 Google LLC
# Licensed under the GNU General Public License v3 (GPLv3).

package WoLFPSort::FastaCleaner;

use strict;
use warnings;
use Carp qw(croak);

our $VERSION = '0.200';

my $HEAD_RE = qr/^\s*>/ax;

sub new {
    my ( $class, %args ) = @_;
    return bless {
        min_len  => defined $args{min_len}  ? $args{min_len}  : 30,
        max_len  => defined $args{max_len}  ? $args{max_len}  : 10000,
        line_len => defined $args{line_len} ? $args{line_len} : 60,
        html     => $args{html}             ? 1               : 0,
    }, $class;
}

sub _extract_id {
    my ( $self, $id_line, $unlabeled_ref ) = @_;
    if ( $id_line =~ /$HEAD_RE \s* (\S+) \s*/ax ) {
        return ( $1, $id_line );
    }
    if ( $id_line =~ /$HEAD_RE \s* $/ax ) {
        ${$unlabeled_ref}++;
        my $id = "queryProtein${$unlabeled_ref}";
        return ( $id, ">$id" );
    }
    croak "Invalid FASTA header line: '$id_line'";
}

sub _flush_record {
    my ( $self, $out_fh, $id, $id_line, $seq ) = @_;
    my $endl = $self->{html} ? "<BR>\n" : "\n";
    my $len  = length $seq;

    if ( $len < $self->{min_len} ) {
        print {$out_fh}
            "#Warning: Sequence $id skipped because length $len less than $self->{min_len}$endl";
        return;
    }
    if ( $len > $self->{max_len} ) {
        print {$out_fh}
            "#Warning: Sequence $id skipped because length $len greater than $self->{max_len}$endl";
        return;
    }
    if ( $seq =~ /O/ax ) {
        print {$out_fh}
            "#Warning: Sequence $id skipped because it contains a non-standard amino acid 'O'$endl";
        print {$out_fh} "#seq is: $seq$endl";
        return;
    }
    if ( $seq !~ /^M/ax ) {
        print {$out_fh} "#Warning: sequence $id does not start with M$endl";
    }

    print {$out_fh} "$id_line\n";
    my $w = $self->{line_len};
    for ( my $pos = 0; $pos < $len; $pos += $w ) {
        print {$out_fh} substr( $seq, $pos, $w ), "\n";
    }
    return;
}

sub _clean_seq_line {
    my ($line) = @_;
    $line =~ tr/a-z/A-Z/;
    $line =~ tr/A-Z//cd;
    return $line;
}

sub clean_fh {
    my ( $self, $in_fh, $out_fh ) = @_;
    croak 'Expected input and output filehandles' unless $in_fh && $out_fh;

    my $endl      = $self->{html} ? "<BR>\n" : "\n";
    my $unlabeled = 0;
    my $seq       = '';
    my $id_line   = '';

    my $line = <$in_fh>;
    while ( defined $line ) {
        $line =~ tr/\015//d;
        last if $line !~ /^\#/ax;
        $line = <$in_fh>;
    }
    croak "Error: empty input$endl" unless defined $line;

    if ( $line =~ /$HEAD_RE/ax ) {
        chomp $line;
        $line =~ s/$HEAD_RE/>/ax;
        $id_line = $line;
    }
    else {
        $id_line = '> queryProtein';
        $seq .= _clean_seq_line($line);
    }
    my $id;
    ( $id, $id_line ) = $self->_extract_id( $id_line, \$unlabeled );

    while ( my $next = <$in_fh> ) {
        $next =~ tr/\015//d;
        next if $next =~ /^\#/ax;
        if ( $next =~ /$HEAD_RE/ax ) {
            $self->_flush_record( $out_fh, $id, $id_line, $seq );
            $seq = '';
            chomp $next;
            $next =~ s/$HEAD_RE/>/ax;
            ( $id, $id_line ) = $self->_extract_id( $next, \$unlabeled );
        }
        else {
            $seq .= _clean_seq_line($next);
        }
    }
    $self->_flush_record( $out_fh, $id, $id_line, $seq );
    return 1;
}

1;

__END__

=head1 NAME

WoLFPSort::FastaCleaner - FASTA input validation and normalization for WoLF PSORT

=head1 SYNOPSIS

    use WoLFPSort::FastaCleaner;

    my $cleaner = WoLFPSort::FastaCleaner->new();
    $cleaner->clean_fh(\*STDIN, \*STDOUT);

=cut
