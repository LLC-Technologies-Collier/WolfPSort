# Copyright 2026 Google LLC
# Licensed under the GNU General Public License v3 (GPLv3).

package WoLFPSort::Backend;

use strict;
use warnings;
use Scalar::Util qw(blessed);
use WoLFPSort;

our $VERSION = '0.200';

sub for_input {
    my ( $class, $input ) = @_;
    if ( blessed($input) && $input->isa('PDL') ) {
        require WoLFPSort::Backend::PDL;
        return 'WoLFPSort::Backend::PDL';
    }
    require WoLFPSort::Backend::PurePerl;
    return 'WoLFPSort::Backend::PurePerl';
}

sub current {
    if ( WoLFPSort::has_pdl() ) {
        require WoLFPSort::Backend::PDL;
        return 'WoLFPSort::Backend::PDL';
    }
    require WoLFPSort::Backend::PurePerl;
    return 'WoLFPSort::Backend::PurePerl';
}

1;

__END__

=head1 NAME

WoLFPSort::Backend - Runtime backend selector (PDL or PurePerl)

=head1 SYNOPSIS

    use WoLFPSort::Backend;

    my $backend = WoLFPSort::Backend->current();

=cut
