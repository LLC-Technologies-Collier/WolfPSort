# Copyright 2026 Google LLC
# Licensed under the GNU General Public License v3 (GPLv3).

use strict;
use warnings;
use ExtUtils::MakeMaker;

WriteMakefile(
    NAME             => 'WoLFPSort',
    AUTHOR           => 'Paul Horton, C.J. Collier',
    VERSION_FROM     => 'lib/WoLFPSort.pm',
    ABSTRACT         => 'Protein Subcellular Localization Prediction',
    LICENSE          => 'gpl_3',
    MIN_PERL_VERSION => '5.016',
    XSMULTI          => 1,
    INC              => '-Isrc',
    LIBS             => ['-lm'],
    PREREQ_PM        => {
        'Carp'         => 0,
        'Getopt::Long' => 0,
        'Scalar::Util' => 0,
    },
    TEST_REQUIRES => {
        'Test::More' => '0.98',
        'File::Temp' => 0,
    },
    EXE_FILES  => [ 'bin/runWolfPsortSummary', 'bin/wolfPredict.pl', ],
    META_MERGE => {
        'meta-spec' => { version => 2 },
        prereqs     => { runtime => { recommends => { 'PDL' => '2.080', }, }, },
        resources   => {
            repository => {
                type => 'git',
                url  => 'https://github.com/LLC-Technologies-Collier/WolfPSort.git',
                web  => 'https://github.com/LLC-Technologies-Collier/WolfPSort',
            },
            bugtracker =>
                { web => 'https://github.com/LLC-Technologies-Collier/WolfPSort/issues', },
        },
    },
);
