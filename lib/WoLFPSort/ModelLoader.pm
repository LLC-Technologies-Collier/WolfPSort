# Copyright (C) 2004-2006 Paul B. Horton
# Copyright 2026 Google LLC
# Licensed under the GNU General Public License v3 (GPLv3).

package WoLFPSort::ModelLoader;

use strict;
use warnings;
use Carp qw(croak);

our $VERSION = '0.200';

sub _slurp_lines {
    my ( $path, $label ) = @_;
    open my $fh, '<', $path or croak "Cannot open $label file '$path': $!";
    my @lines = <$fh>;
    close $fh or croak "Cannot close $label file '$path': $!";
    return \@lines;
}

sub load_weights {
    my ($path) = @_;
    my $lines = _slurp_lines( $path, 'weight' );

    my ( @features, @w1, @w2, $k );
    my %handlers = (
        feature => sub { @features = split /\s+/ax, $_[0] },
        w1      => sub { @w1       = split /\s+/ax, $_[0] },
        w2      => sub { @w2       = split /\s+/ax, $_[0] },
    );

    for my $line ( @{$lines} ) {
        $line =~ s/\#.*//ax;
        $line =~ s/^\s+|\s+$//agx;
        next unless length $line;

        if ( $line =~ /^(\w+)\s+(.*)$/ax && exists $handlers{$1} ) {
            $handlers{$1}->($2);
        }
        elsif ( $line =~ /^k\s*=\s*(\d+)/ax ) {
            $k = int $1;
        }
    }

    croak "Weight file '$path' missing feature/w1/w2/k definitions"
        unless @features && @w1 && @w2 && defined $k;
    croak "Weight vector length mismatch in '$path'"
        unless @w1 == @features && @w2 == @features;

    return {
        feature_names => \@features,
        w1            => [ map { $_ + 0 } @w1 ],
        w2            => [ map { $_ + 0 } @w2 ],
        k             => $k,
    };
}

sub load_utility {
    my ($path) = @_;
    my $lines = _slurp_lines( $path, 'utility' );

    my ( @classes, @rows );
    for my $line ( @{$lines} ) {
        $line =~ s/\#.*//ax;
        $line =~ s/^\s+|\s+$//agx;
        next unless length $line;

        if ( !@classes ) {
            @classes = split /\s+/ax, $line;
        }
        else {
            my @vals = map { $_ + 0 } split /\s+/ax, $line;
            push @rows, \@vals;
        }
    }

    my $n_cls = scalar @classes;
    if ( $n_cls <= 0 || @rows != $n_cls ) {
        croak "Utility matrix in '$path' is empty or non-square";
    }

    my %class_idx;
    for my $i ( 0 .. $n_cls - 1 ) {
        $class_idx{ $classes[$i] } = $i;
    }

    return {
        class_names => \@classes,
        class_index => \%class_idx,
        matrix      => \@rows,
    };
}

sub load_training {
    my ( $path, $feature_names, $class_index, %opts ) = @_;
    my $pack_doubles = $opts{pack_doubles} ? 1 : 0;
    my $lines        = _slurp_lines( $path, 'training' );

    my $header = shift @{$lines};
    croak "Training file '$path' is empty" unless defined $header;
    $header =~ s/\#.*//ax;
    $header =~ s/^\s+|\s+$//agx;
    my @hdr = split /\s+/ax, $header;
    splice @hdr, 0, 2;

    croak "Training feature count mismatch in '$path'"
        unless @hdr == @{$feature_names};

    my ( @ids, @classes, @class_ids, @rows );
    my $buf = '';

    for my $line ( @{$lines} ) {
        $line =~ s/\#.*//ax;
        $line =~ s/^\s+|\s+$//agx;
        next unless length $line;

        my @fields = split /\s+/ax, $line;
        my $id     = shift @fields;
        my $cls    = shift @fields;
        croak "Unknown class '$cls' for training instance '$id'"
            unless exists $class_index->{$cls};

        push @ids,       $id;
        push @classes,   $cls;
        push @class_ids, $class_index->{$cls};
        if ($pack_doubles) {
            $buf .= pack 'd*', @fields;
        }
        else {
            push @rows, [ map { $_ + 0 } @fields ];
        }
    }

    my $n_rows = scalar @ids;
    croak "No training instances found in '$path'" unless $n_rows;

    return {
        ids       => \@ids,
        classes   => \@classes,
        class_ids => \@class_ids,
        n_rows    => $n_rows,
        n_feat    => scalar( @{$feature_names} ),
        ( $pack_doubles ? ( packed_buf => $buf ) : ( rows => \@rows ) ),
    };
}

1;

__END__

=head1 NAME

WoLFPSort::ModelLoader - Parsers for WoLF PSORT weight, utility, and training files

=head1 SYNOPSIS

    use WoLFPSort::ModelLoader;

    my $w = WoLFPSort::ModelLoader::load_weights('data/fungi.wolfw');
    my $u = WoLFPSort::ModelLoader::load_utility('data/fungi.wolfu');
    my $t = WoLFPSort::ModelLoader::load_training(
        'data/fungi.wolff', $w->{feature_names}, $u->{class_index},
    );

=cut
