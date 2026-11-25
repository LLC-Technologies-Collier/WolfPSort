# Copyright (C) 2004-2006 Paul B. Horton
# Copyright 2026 Google LLC
# Licensed under the GNU General Public License v3 (GPLv3).

package WoLFPSort::CLI::Predict;

use strict;
use warnings;
use Carp         qw(croak);
use Getopt::Long qw(GetOptionsFromArray :config no_ignore_case);
use WoLFPSort::KNN;

our $VERSION = '0.200';

my $USAGE = <<'USAGE';
Usage: wolfPredict.pl [OPTIONS] weightFile trainingDataFile utilityMatrixFile

  OPTIONS

    --ignore-self
        If instance in training data has the same id as the query id
        ignore that instance. Useful for jackknife accuracy testing.

    --print-instance-comment
        Print the instance comment (and a newline) instead of the
        instance id (and space character) before the prediction result.

    --print-neighbors
        Print the nearest neighbors used to make each prediction.

    --print-all-scores
        Print localization class and utility for all classes.
USAGE

sub filter_scores {
    my ( $scores, $print_all ) = @_;
    return @{$scores} if $print_all;

    my @kept;
    my $top       = @{$scores} ? $scores->[0]{score} : 0;
    my $min_score = $top >= 10 ? 2.5                 : 1.0;

    for my $entry ( @{$scores} ) {
        last if $entry->{score} < $min_score;
        push @kept, $entry;
    }
    return @kept ? @kept : ( @{$scores} ? ( $scores->[0] ) : () );
}

sub format_number {
    my ($val) = @_;
    return '0' if abs($val) < 1e-9;
    my $s = sprintf '%.3f', $val;
    $s =~ s/\.?0+$//ax;
    return $s;
}

sub predict_stream {
    my ( $knn, $in_fh, $out_fh, %opts ) = @_;
    my $ignore_self            = $opts{ignore_self}            ? 1 : 0;
    my $print_instance_comment = $opts{print_instance_comment} ? 1 : 0;
    my $print_neighbors        = $opts{print_neighbors}        ? 1 : 0;
    my $print_all_scores       = $opts{print_all_scores}       ? 1 : 0;

    print {$out_fh} '# k used for kNN is: ', $knn->k, "\n";

    my $header = <$in_fh>;
    return 0 unless defined $header;

    while ( my $line = <$in_fh> ) {
        my $comment = '';
        if ( $line =~ s/\#(.*)$//ax ) {
            $comment = $1;
            $comment =~ s/\r?\n$//ax;
        }
        $line =~ s/^\s+|\s+$//agx;
        next unless length $line;

        my @fields = split /\s+/ax, $line;
        my $id     = shift @fields;
        my $cls    = shift @fields;

        my $res = $knn->predict_vector(
            \@fields,
            id          => $id,
            ignore_self => $ignore_self,
        );

        my @display  = filter_scores( $res->{scores}, $print_all_scores );
        my $pred_str = join ', ',
            map { $_->{class} . q{ } . format_number( $_->{score} ) } @display;

        if ($print_instance_comment) {
            print {$out_fh} "$comment\n$pred_str\n";
        }
        else {
            print {$out_fh} "$id $pred_str\n";
        }

        if ($print_neighbors) {
            for my $nb ( @{ $res->{neighbors} } ) {
                printf {$out_fh} "  %-11s %-10s %s\n",
                    $nb->{id},
                    $nb->{class},
                    format_number( $nb->{similarity} );
            }
        }
    }

    return 0;
}

sub run {
    my ( $argv, $in_fh, $out_fh ) = @_;
    $in_fh  ||= \*STDIN;
    $out_fh ||= \*STDOUT;

    my @args = @{$argv};
    my %opts = (
        ignore_self            => 0,
        print_instance_comment => 0,
        print_neighbors        => 0,
        print_all_scores       => 0,
    );

    GetOptionsFromArray(
        \@args,
        'ignore-self'            => \$opts{ignore_self},
        'print-instance-comment' => \$opts{print_instance_comment},
        'print-neighbors'        => \$opts{print_neighbors},
        'print-all-scores'       => \$opts{print_all_scores},
    ) or croak $USAGE;

    @args == 3 or croak $USAGE;
    my ( $weight_file, $training_file, $utility_file ) = @args;

    my $knn = WoLFPSort::KNN->new(
        weight_file   => $weight_file,
        training_file => $training_file,
        utility_file  => $utility_file,
    );

    return predict_stream( $knn, $in_fh, $out_fh, %opts );
}

1;

__END__

=head1 NAME

WoLFPSort::CLI::Predict - CLI driver for WoLF PSORT kNN prediction

=head1 SYNOPSIS

    use WoLFPSort::CLI::Predict;
    exit WoLFPSort::CLI::Predict::run(\@ARGV, \*STDIN, \*STDOUT);

=cut
