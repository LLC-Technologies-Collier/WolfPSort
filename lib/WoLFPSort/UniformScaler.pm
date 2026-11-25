# Copyright (C) 2004-2006 Paul B. Horton
# Copyright 2026 Google LLC
# Licensed under the GNU General Public License v3 (GPLv3).

package WoLFPSort::UniformScaler;

use strict;
use warnings;
use Carp qw(croak);
use WoLFPSort::Backend;

our $VERSION = '0.200';

sub new {
    my ( $class, $matrix ) = @_;
    my $backend = WoLFPSort::Backend->for_input($matrix);
    my ( $tables, $n_features, $n_samples ) = $backend->build_scaler_tables($matrix);

    return bless {
        backend    => $backend,
        tables     => $tables,
        n_features => $n_features,
        n_samples  => $n_samples,
    }, $class;
}

sub n_features {
    my ($self) = @_;
    return $self->{n_features};
}

sub n_samples {
    my ($self) = @_;
    return $self->{n_samples};
}

sub transform {
    my ( $self, $X ) = @_;
    return $self->{backend}->transform_scaler( $self->{tables}, $self->{n_features}, $X );
}

1;

__END__

=head1 NAME

WoLFPSort::UniformScaler - Empirical-CDF feature scaler for WoLF PSORT

=head1 SYNOPSIS

    use WoLFPSort::UniformScaler;

    my $scaler = WoLFPSort::UniformScaler->new($train_matrix);
    my $scaled_query = $scaler->transform($query_vector);

=head1 DESCRIPTION

Maps each feature column to uniform C<[0, 100]> mid-rank percentiles with
out-of-sample linear interpolation/extrapolation (Horton et al., APBC06,
Section 2.3.1, p. 3), delegating numerical operations to L<WoLFPSort::Backend>.

=cut
