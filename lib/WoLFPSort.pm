# Copyright (C) 2004-2006 Paul B. Horton
# Copyright 2026 Google LLC
# Licensed under the GNU General Public License v3 (GPLv3).

package WoLFPSort;

use strict;
use warnings;

our $VERSION = '0.200';

my $_HAS_PDL;
my $_HAS_XS;

sub has_pdl {
    return 0 if $ENV{WOLFPSORT_NO_PDL};
    if ( !defined $_HAS_PDL ) {
        $_HAS_PDL = eval { require PDL::Lite; 1 } ? 1 : 0;
    }
    return $_HAS_PDL;
}

sub has_xs {
    return 0 if $ENV{WOLFPSORT_NO_XS};
    if ( !defined $_HAS_XS ) {
        $_HAS_XS = eval { require WoLFPSort::XS; 1 } ? 1 : 0;
    }
    return $_HAS_XS;
}

1;

__END__

=head1 NAME

WoLFPSort - Protein Subcellular Localization Prediction

=head1 SYNOPSIS

    use WoLFPSort;
    use WoLFPSort::KNN;

    my $knn = WoLFPSort::KNN->new(
        weight_file   => 'data/fungi.wolfw',
        training_file => 'data/fungi.wolff',
        utility_file  => 'data/fungi.wolfu',
    );

=head1 DESCRIPTION

C<WoLFPSort> is a modular Perl implementation of Paul Horton's WoLF PSORT
subcellular localization predictor. When L<PDL> is installed, numerical
routines use PDL ndarrays; otherwise (or when C<WOLFPSORT_NO_PDL=1> is set)
they run in pure Perl.

=head1 MODULES

=over 4

=item * L<WoLFPSort::UniformScaler> - Empirical-CDF quantile feature scaling

=item * L<WoLFPSort::KNN> - Weighted hybrid L1/L2 I<k>NN classifier with utility scoring

=item * L<WoLFPSort::AminoAcidIndex> - Amino acid index and sliding-window sums

=item * L<WoLFPSort::AminoAcidIndexer> - Amino acid code-to-index mapping and sequence encoding

=item * L<WoLFPSort::FastaCleaner> - FASTA input validation and normalization

=item * L<WoLFPSort::MotifMatch>, L<WoLFPSort::RegexMotif>, L<WoLFPSort::RegexMotifs> - Regular-expression motif matching

=item * L<WoLFPSort::PSORTII::KNN> - Legacy PSORT II Euclidean I<k>NN

=item * L<WoLFPSort::CLI::Predict>, L<WoLFPSort::CLI::Summary> - CLI drivers

=back

=head1 TYPES

When L<PDL> is enabled, all numerical computations use real-valued C<double>
L<PDL> ndarrays (with C<long> index vectors). Complex PDL types are not
supported because protein sequence indices, empirical quantiles, and
localization utility scores are strictly real-valued.

=head1 LICENSE

GNU General Public License v3 (GPLv3).

=cut
