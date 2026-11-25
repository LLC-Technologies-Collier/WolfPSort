# Copyright (C) 2004-2006 Paul B. Horton
# Copyright 2026 Google LLC
# Licensed under the GNU General Public License v3 (GPLv3).

package WoLFPSort::XS;

use strict;
use warnings;
use XSLoader;

our $VERSION = '0.200';

XSLoader::load( 'WoLFPSort::XS', $VERSION );

1;

__END__

=head1 NAME

WoLFPSort::XS - In-process XS bindings to the C99 WoLF PSORT kNN core

=head1 SYNOPSIS

    use WoLFPSort::XS;

    my $model = WoLFPSort::XS->new(
        'data/fungi.wolfw',
        'data/fungi.wolff',
        'data/fungi.wolfu',
    );
    my $res = $model->predict_vector(\@features, 'RCC1_YEAST', 0);

=cut
