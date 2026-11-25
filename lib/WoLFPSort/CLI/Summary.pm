package WoLFPSort::CLI::Summary;
#
#  Author: Paul Horton (remnants of code from CJ Collier may remain).
#  Maintainer: Paul Horton.
#  Organization: Computational Biology Research Center, AIST, Japan
#  Creation Date: Summer 2004
#  Last Modification $Date: 2006/08/31 07:21:44 $
#  Copyright: Paul Horton
#  Copyright 2026 Google LLC
#  Licensing: Anyone may distribute this under the GNU public license. For other
#             license information contact: horton-p AT aist.go.jp
#
#  Description: runWolfPsortSummary runs the programs required for
#  predicting localization sites of fasta queries. Only a
#  simple summary of the predictions is output. Other scripts are available
#  for more elaborate html output.
#
#  Input: (multi)fasta query
#
#  Output: Localization site predictions
#
#  Args: organism type from {fungi, plant, animal}
#
#  Dependencies: This program depends on some files included in the
#                installation in its directory or in ../data.
#

use strict;
use warnings;
use Carp qw(croak);
use File::Spec;
use Getopt::Long qw(GetOptionsFromArray :config posix_default bundling permute);
use Pod::Usage   qw(pod2usage);
use WoLFPSort::FastaCleaner;
use WoLFPSort::FeatureExtractor;
use WoLFPSort::KNN;
use WoLFPSort::CLI::Predict;

our $VERSION = '0.200';

sub _resolve_pod_path {
    my ( $bin_dir, $pod_path ) = @_;
    if ( defined $ENV{LANG} && $ENV{LANG} =~ /^ja_JP/ax ) {
        my $ja_pod
            = File::Spec->catfile( $bin_dir, 'pod', 'ja_JP', 'runWolfPsortSummary.ja_JP.eucjp.pod',
            );
        return $ja_pod if -e $ja_pod;
    }
    return $pod_path;
}

sub _parse_args {
    my ( $argv, $pod_path, $out_fh ) = @_;
    my @args = @{$argv};
    my %opts = (
        just_print            => 0,
        just_compute_features => 0,
        append_head           => 0,
        print_neighbors       => 0,
    );
    my ( $man, $help, $usage ) = ( 0, 0, 0 );

    my $ok = GetOptionsFromArray(
        \@args,
        'n|just-print'                  => \$opts{just_print},
        'f|just-compute-features'       => \$opts{just_compute_features},
        'd|print-sequence-descriptions' => \$opts{append_head},
        'print-neighbors'               => \$opts{print_neighbors},
        'man'                           => \$man,
        'help'                          => \$help,
        'usage'                         => \$usage,
    );
    $ok || pod2usage( -verbose => 0, -input => $pod_path );
    $man   && pod2usage( -verbose => 2, -input => $pod_path );
    $help  && pod2usage( -verbose => 1, -input => $pod_path );
    $usage && pod2usage( -verbose => 0, -input => $pod_path );

    my $organism = shift @args;
    if ( !defined $organism || $organism !~ /^(?:fungi|plant|animal)$/ax ) {
        print {$out_fh} "OrganismType should be one of {\"fungi\", \"plant\", \"animal\"}\n"
            . "  pipe fasta query sequences in from standard in\n\n";
        pod2usage( -verbose => 0, -input => $pod_path );
    }
    if (@args) {
        print {$out_fh}
            "\nCommand Line Parsing Error; Do not know what to do with argument \"$args[0]\"\n\n";
        pod2usage( -verbose => 0, -input => $pod_path );
    }

    $opts{organism} = $organism;
    return \%opts;
}

sub run {
    my ( $argv, %ctx ) = @_;
    my $bin_dir  = $ctx{bin_dir} || croak 'Missing bin_dir';
    my $pod_path = _resolve_pod_path( $bin_dir, $ctx{pod_path} || $0 );
    my $in_fh    = $ctx{in_fh}  || \*STDIN;
    my $out_fh   = $ctx{out_fh} || \*STDOUT;

    my $opts      = _parse_args( $argv, $pod_path, $out_fh );
    my $organism  = $opts->{organism};
    my $data_dir  = File::Spec->catdir( $bin_dir, File::Spec->updir, 'data' );
    my $psort_dir = File::Spec->catdir( $bin_dir, 'psortModifiedForWolfFiles' );

    if ( $opts->{just_print} ) {
        print {$out_fh}
            "WoLFPSort::FastaCleaner -> WoLFPSort::FeatureExtractor -> WoLFPSort::KNN ($organism)\n";
        return 0;
    }

    my $cleaned = '';
    {
        open my $clean_fh, '>', \$cleaned or croak "Cannot open scalar buffer: $!";
        my $cleaner = WoLFPSort::FastaCleaner->new();
        $cleaner->clean_fh( $in_fh, $clean_fh );
        close $clean_fh or croak "Cannot close scalar buffer: $!";
    }

    my $features = '';
    {
        open my $cfh, '<', \$cleaned  or croak "Cannot open cleaned buffer: $!";
        open my $ffh, '>', \$features or croak "Cannot open feature buffer: $!";
        my $fe = WoLFPSort::FeatureExtractor->new(
            psort_dir   => $psort_dir,
            append_head => $opts->{append_head},
        );
        $fe->extract_stream( $cfh, $ffh );
        close $cfh or croak "Cannot close cleaned buffer: $!";
        close $ffh or croak "Cannot close feature buffer: $!";
    }

    if ( $opts->{just_compute_features} ) {
        print {$out_fh} $features;
        return 0;
    }

    my $knn = WoLFPSort::KNN->new(
        weight_file   => File::Spec->catfile( $data_dir, "$organism.wolfw" ),
        training_file => File::Spec->catfile( $data_dir, "$organism.wolff" ),
        utility_file  => File::Spec->catfile( $data_dir, "$organism.wolfu" ),
    );

    open my $feat_fh, '<', \$features or croak "Cannot open feature buffer: $!";
    my $rc = WoLFPSort::CLI::Predict::predict_stream(
        $knn,
        $feat_fh,
        $out_fh,
        print_instance_comment => $opts->{append_head},
        print_neighbors        => $opts->{print_neighbors},
    );
    close $feat_fh or croak "Cannot close feature buffer: $!";
    return $rc;
}

1;

__END__

=head1 NAME

WoLFPSort::CLI::Summary - In-process CLI driver for runWolfPsortSummary

=head1 SYNOPSIS

    use WoLFPSort::CLI::Summary;
    exit WoLFPSort::CLI::Summary::run(\@ARGV, bin_dir => $bin_dir, pod_path => $0);

=cut
