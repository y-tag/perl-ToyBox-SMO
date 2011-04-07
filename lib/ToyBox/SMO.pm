package ToyBox::SMO;

use strict;
use warnings;

use List::Util;

our $VERSION = '0.0.1';

sub new {
    my $class = shift;
    my $self = {
        dnum => 0,
        data => [],
        b => 0.0,
        C => 1.0,
    } ;
    bless $self, $class;
}

sub add_instance {
    my ($self, %params) = @_;

    my $attributes = $params{attributes} or die "No params: attributes";
    my $label      = $params{label}     or die "No params: label";
    die "attributes is not hash ref"   unless ref($attributes) eq 'HASH';
    die "attributes is empty hash ref" unless keys %$attributes;

    my %copy_attr = %$attributes;

    my $i = $self->{dnum};

    my $datum = {i => $i,
                 x => \%copy_attr,
                 y => $label,
                 alpha => 0};

    push(@{$self->{data}}, $datum);
    $self->{dnum}++;
}

sub train {
    my ($self, %params) = @_;

    my $C = $params{C};
    $C = 1.0 unless defined($C);
    die "C is le 0" unless scalar($C) > 0;

    my $progress_cb = $params{progress_cb};
    my $verbose = 0;
    $verbose = 1 if defined($progress_cb) && $progress_cb eq 'verbose';

    my $num_changed = 0;
    my $examine_all = 1;

    $self->{C} = $C;

    foreach my $datum (@{$self->{data}}) {
        $datum->{alpha} = $C / 2;
    }

    while ($num_changed > 0 || $examine_all) {
        $num_changed = 0;
        if ($examine_all) {
            my @data = @{$self->{data}};
            foreach my $datum (@data) {
                $num_changed += $self->examine_example($datum);
            }
        } else {
            my @data = grep {$_->{alpha} > 0 && $_->{alpha} < $C} @{$self->{data}};
            foreach my $datum (@data) {
                $num_changed += $self->examine_example($datum);
            }
        }

        if ($examine_all == 1) {
            $examine_all = 0;
        } elsif ($num_changed == 0) {
            $examine_all = 1;
        }
    }
    1;
}

sub examine_example {
    my ($self, $datum2) = @_;

    my $y2 = $datum2->{y};
    my $alpha2 = $datum2->{alpha};
    my $E2 = $self->calc_error($datum2);
    my $r2 = $E2 * $y2;

    $datum2->{error} = $E2;

    my $C = $self->{C};
    my $tol = 10e-15;

    if (($r2 < -$tol && $alpha2 < $C) || ($r2 > $tol && $alpha2 > 0)) {
        my @data = grep {$_->{alpha} > 0 && $_->{alpha} < $C } @{$self->{data}};
        if (@data > 1) {
            my $datum1 = $self->second_choice($datum2);
            return 1 if $self->take_step($datum1, $datum2);
        }

        foreach my $datum1 (List::Util::shuffle @data) {
           return 1 if  $self->take_step($datum1, $datum2);
        }

        foreach my $datum1 (List::Util::shuffle @{$self->{data}}) {
           return 1 if  $self->take_step($datum1, $datum2);
        }
    }
    return 0;
}

sub take_step {
    my ($self, $datum1, $datum2) = @_;

    return 0 if $datum1->{i} == $datum2->{i};

    my $C = $self->{C};
    my $eps = 10e-5;

    my $alpha1 = $datum1->{alpha};
    my $y1 = $datum1->{y};
    my $E1 = $self->calc_error($datum1);

    my $alpha2 = $datum2->{alpha};
    my $y2 = $datum2->{y};
    my $E2 = $self->calc_error($datum2);

    my $s = $y1 * $y2;

    my ($L, $H);

    if ($y1 == $y2) {
        $L = List::Util::max(0, $alpha2 + $alpha1 - $C);
        $H = List::Util::min($C, $alpha2 + $alpha1);
    } else {
        $L = List::Util::max(0, $alpha2 - $alpha1);
        $H = List::Util::min($C, $C + $alpha2 - $alpha1);
    }

    return 0 if $L == $H;

    my $k11 = $self->kernel($datum1->{x}, $datum1->{x});
    my $k12 = $self->kernel($datum1->{x}, $datum2->{x});
    my $k22 = $self->kernel($datum2->{x}, $datum2->{x});
    my $eta = $k11 + $k22 - 2 * $k12;

    my $a2;

    if ($eta > 0) {
        $a2 = $alpha2 + $y2 * ($E1 - $E2) / $eta;
        if    ($a2 < $L) { $a2 = $L; }
        elsif ($a2 > $H) { $a2 = $H; }
    } else {
        my $f1 = $y1*($E1 + $b) - $alpha1*$k11 - $s*$alpha2*$k12;
        my $f2 = $y2*($E2 + $b) - $s*$alpha1*$k12 - $alpha2*$k22;
        my $L1 = $alpha1 + $s*($alpha2 - $L);
        my $H1 = $alpha1 + $s*($alpha2 - $H);
        my $Lobj = $L1*$f1 + $L*$f2 + $L1*$L1*$k11/2 + $L*$L*$k22/2 + $s*$L* $L1*$k12;
        my $Hobj = $H1*$f1 + $H*$f2 + $H1*$H1*$k11/2 + $H*$H*$k22/2 + $s*$H* $H1*$k12;

        if ($Lobj < $Hobj - $eps) {
            $a2 = $L;
        } elsif ($Lobj > $Hobj + $eps) {
            $a2 = $H;
        } else {
            $a2 = $alpha2;
        }
    }

    if (abs($a2-$alpha2) < $eps*($a2+$alpha2+$eps)) {
        return 0;
    }

    my $a1 = $alpha1 + $s * ($alpha2 - $a2);

    my $b = $self->{b};
    my $b1 = $E1 + $y1*($a1-$alpha1)*$k11 + $y2*($a2-$alpha2)*$k12 + $b;
    my $b2 = $E2 + $y1*($a1-$alpha1)*$k12 + $y2*($a2-$alpha2)*$k22 + $b;
    $self->{b} = ($b1 + $b2) / 2;

    $datum1->{alpha} = $a1;
    $datum2->{alpha} = $a2;

    $self->update_error();

    return 1;
}

sub calc_error {
    my ($self, $datum) = @_;

    my $E = $self->predict(attributes => $datum->{x}) - $datum->{y};
    $datum->{error} = $E;

    $E;
}

sub update_error {
    my $self = shift;

    foreach my $datum (@{$self->{data}}) {
        $self->calc_error($datum);
    }
}

sub second_choice {
    my ($self, $datum2) = @_;

    my $E2 = $datum2->{error};
    my $datum1 = $self->{data}[0];
    my $E1 = defined($datum1->{error}) ? $datum1->{error} : $self->calc_error($datum1);

    if ($E2 > 0) {
        foreach my $datum (@{$self->{data}}) {
            $self->calc_error($datum) unless defined($datum->{error});
            if ($E1 > $datum->{error}) {
                $E1 = $datum->{error};
                $datum1 = $datum;
            }
        }
    } else {
        foreach my $datum (@{$self->{data}}) {
            $self->calc_error($datum) unless defined($datum->{error});
            if ($E1 < $datum->{error}) {
                $E1 = $datum->{error};
                $datum1 = $datum;
            }
        }
    }

    $datum1;
}

sub kernel {
    my ($self, $x1, $x2) = @_;

    my $ret = 0;
    while (my ($f, $v) = each %$x1) {
        $ret += $v * $x2->{$f} if defined($x2->{$f});
    }
    return $ret;
}


sub predict {
    my ($self, %params) = @_;

    my $attributes = $params{attributes} or die "No params: attributes";
    die "attributes is not hash ref"   unless ref($attributes) eq 'HASH';
    die "attributes is empty hash ref" unless keys %$attributes;

    my $score = 0;
    foreach my $datum (@{$self->{data}}) {
        $score += $datum->{y} * $datum->{alpha} * $self->kernel($datum->{x}, $attributes) unless $datum->{alpha} == 0;
    }
    $score -= $self->{b};

    $score;
}

1;
__END__


=head1 NAME

ToyBox::SMO - Support Vector Machine with Sequential Minimal Optimization

=head1 SYNOPSIS

  use ToyBox::SMO;

  my $smo = ToyBox::SMO->new();
  
  $smo->add_instance(
      attributes => {a => 2, b => 3},
      label => 1 
  );
  
  $smo->add_instance(
      attributes => {c => 3, d => 1},
      label => -1
  );
  
  $smo->train();
  
  my $score = $smo->predict(
                  attributes => {a => 1, b => 1, d => 1, e =>1}
              );

=head1 DESCRIPTION

=head1 AUTHOR

TAGAMI Yukihiro <tagami.yukihiro@gmail.com>

=head1 LICENSE

This library is distributed under the term of the MIT license.

L<http://opensource.org/licenses/mit-license.php>

