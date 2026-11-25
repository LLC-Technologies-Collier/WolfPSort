# Copyright (C) 2004-2006 Paul B. Horton
# Copyright 2026 Google LLC
# Licensed under the GNU General Public License v3 (GPLv3).

package WoLFPSort::FeatureExtractor;
## no critic (Variables::ProhibitPackageVars)

use strict;
use warnings;
use Carp qw(croak);
use File::Spec;
use WoLFPSort::AminoAcidIndex;
use WoLFPSort::RegexMotifs;

our $VERSION = '0.200';

my @AMINO_ACID_CHARS = qw( A R N D C Q E G H I L K M F P S T W Y V );

my @VARIABLE_NAMES = sort qw(
    act alm bac caa dna erl erm gpi gvh leu
    m1a m1b m2 m3a m3b mNt m_ mip mit myr nuc
    pox psg px2 rib rnp tms top tyr vac yqr
    maxHydropathy5_29_11 maxHydropathy0_29_12
    averageNegativeCharge0_24 maxNegativeCharge0_19_12
);

my $_LEGACY_LOADED = 0;

sub new {
    my ( $class, %args ) = @_;
    my $psort_dir = $args{psort_dir} || croak 'Missing required argument psort_dir';
    my $self      = bless {
        psort_dir   => $psort_dir,
        motif_dir   => File::Spec->catdir( $psort_dir, 'motifLibraries' ),
        min_seq_len => defined $args{min_seq_len} ? $args{min_seq_len} : 30,
        append_head => $args{append_head}         ? 1                  : 0,
        motif_cache => {},
    }, $class;

    $self->_load_legacy_psort();
    return $self;
}

sub _load_legacy_psort {
    my ($self) = @_;
    return if $_LEGACY_LOADED;

    my $mod_dir = File::Spec->catdir( $self->{psort_dir}, 'psortModules' );
    local @INC = ( $mod_dir, @INC );
    require PsortSeq;
    require File::Spec->catfile( $self->{psort_dir}, 'sub1.pl' );
    require File::Spec->catfile( $self->{psort_dir}, 'sub2.pl' );
    $_LEGACY_LOADED = 1;
    return;
}

sub _load_motif_library {
    my ( $self, $name ) = @_;
    if ( !exists $self->{motif_cache}{$name} ) {
        my $path = File::Spec->catfile( $self->{motif_dir}, "$name.regexMotifs" );
        open my $fh, '<', $path or croak "Cannot open motif library '$path': $!";
        $self->{motif_cache}{$name} = WoLFPSort::RegexMotifs->new($fh);
        close $fh or croak "Cannot close motif library '$path': $!";
    }
    return $self->{motif_cache}{$name};
}

# Horton et al. (APBC06), Section 2.3.1 (Candidate Features), p. 3, para. 1 &
# Table 5, p. 7: extracts PSORT/PSORT II sorting signals (gvh, psg, mip, mit,
# alm, m1a/m1b/m3a/m3b/mNt, tms, erl, leu, nuc, pox, tyr, yqr, vac, dna, rib,
# myr, rnp, act), iPSORT hydropathy/charge windows, amino acid composition, and
# sequence length.
sub _compute_seqana {
    my ($self) = @_;
    my %vars;

    $vars{psg} = sub1::psg();
    my ( $gvh, $gvhpos ) = sub1::gvh('eukaryote');
    $vars{gvh} = $gvh;

    my ( $sig, $npos, $middle );
    if ( $vars{psg} > 0 && $vars{gvh} > -2.1 ) {
        $sig = 'cleavable';
        ( $vars{tms}, $vars{alm}, $npos ) = sub1::alom2( $gvhpos + 1, 0.5, -2.0 );
        $middle = int( $gvhpos / 2 );
    }
    else {
        $sig = 'nosignal';
        ( $vars{tms}, $vars{alm}, $npos ) = sub1::alom2( 1, 0.5, -2.0 );
        $middle = $npos + 7;
    }

    $vars{top} = ( $vars{psg} <= 0 && $vars{tms} == 0 ) ? 0 : sub1::mtop($middle);
    my ( $mtype, $start, $end ) = sub1::mtype_assign( $vars{tms}, $sig, $vars{top} );
    croak "$PsortSeq::id: inconsistent ($mtype) $start : $end" if $start > $end;

    $vars{gpi} = ( $vars{tms} > 0 && $mtype eq '1a' && ( $end - $start < 10 ) ) ? 1 : 0;

    $vars{mit} = sub2::mitdisc();
    my ($mip) = sub2::gavel();
    $vars{mip} = $mip;
    $vars{nuc} = sub2::nucdisc();
    $vars{erl} = sub2::hdel();
    $vars{erm} = sub2::erm($mtype);
    $vars{pox} = sub2::pts1();
    $vars{px2} = sub2::pts2();
    $vars{vac} = sub2::vaccalc();
    $vars{rnp} = sub2::rnp1();
    $vars{act} = sub2::actin();
    $vars{myr} = sub2::nmyr();
    $vars{caa} = sub2::isoprenyl();
    $vars{yqr} = sub2::yqrl( $vars{tms}, $start );
    $vars{tyr} = sub2::tyros( $vars{tms}, $start, $end );
    $vars{leu} = sub2::dileu( $vars{tms}, $start, $end );

    $vars{dna} = $self->_load_motif_library('dnaBinding')->num_matches($PsortSeq::seq);
    $vars{rib} = $self->_load_motif_library('ribosomal')->num_matches($PsortSeq::seq);
    $vars{bac} = $self->_load_motif_library('bacterialDnaBinding')->num_matches($PsortSeq::seq);

    my %mtype_map = (
        '1a' => 'm1a',
        '1b' => 'm1b',
        '2 ' => 'm2',
        'Nt' => 'mNt',
        '3a' => 'm3a',
        '3b' => 'm3b',
        '__' => 'm_',
    );
    @vars{qw(m1a m1b m2 mNt m3a m3b m_)} = (0) x 7;
    croak "Undefined mtype value: '$mtype'" unless exists $mtype_map{$mtype};
    $vars{ $mtype_map{$mtype} } = 1;

    my $kd    = WoLFPSort::AminoAcidIndex->new('KYTJ820101');
    my $klein = WoLFPSort::AminoAcidIndex->new('KLEP840101');
    $vars{maxHydropathy5_29_11}      = $kd->maxWindowAverage( \@PsortSeq::seq, 11, 5, 24 );
    $vars{maxHydropathy0_29_12}      = $kd->maxWindowAverage( \@PsortSeq::seq, 12, 0, 29 );
    $vars{averageNegativeCharge0_24} = $klein->seqAverage( \@PsortSeq::seq, 0, 24 );
    $vars{maxNegativeCharge0_19_12}  = $klein->maxWindowAverage( \@PsortSeq::seq, 7, 0, 19 );

    return \%vars;
}

sub extract_stream {
    my ( $self, $in_fh, $out_fh ) = @_;
    croak 'Expected input and output filehandles' unless $in_fh && $out_fh;

    local $main::OPT_V = 0;
    local $main::endl  = "\n";
    $PsortSeq::fastaReader = PsortFastaReader->new();

    print {$out_fh} 'id class ';
    for my $aa (@AMINO_ACID_CHARS) {
        print {$out_fh} " ${aa}cont";
    }
    print {$out_fh} ' length ';
    for my $var (@VARIABLE_NAMES) {
        print {$out_fh} "$var ";
    }
    print {$out_fh} "\n";

    while ( PsortSeq::readOneSequence($in_fh) ) {
        next if $PsortSeq::len < $self->{min_seq_len};
        next if $PsortSeq::badAminoAcidCharsFoundList;

        my $vars    = $self->_compute_seqana();
        my $num_std = $PsortSeq::seq =~ tr/ARNDCQEGHILKMFPSTWYV/ARNDCQEGHILKMFPSTWYV/;

        print {$out_fh} "$PsortSeq::id unknown";
        my %aa_counts;
        for my $ch (@PsortSeq::seq) {
            $aa_counts{$ch}++;
        }
        for my $aa (@AMINO_ACID_CHARS) {
            my $frac = ( $aa_counts{$aa} || 0 ) / $num_std;
            print {$out_fh} " $frac";
        }
        print {$out_fh} " $PsortSeq::len";
        for my $var (@VARIABLE_NAMES) {
            print {$out_fh} " $vars->{$var}";
        }
        print {$out_fh} "#$PsortSeq::fastaHead" if $self->{append_head};
        print {$out_fh} "\n";
    }
    return 1;
}

1;

__END__

=head1 NAME

WoLFPSort::FeatureExtractor - In-process PSORT/WoLF protein sequence feature extractor

=head1 SYNOPSIS

    use WoLFPSort::FeatureExtractor;

    my $fe = WoLFPSort::FeatureExtractor->new(
        psort_dir => 'bin/psortModifiedForWolfFiles',
    );
    $fe->extract_stream($in_fh, $out_fh);

=cut
