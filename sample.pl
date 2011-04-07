#/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;

use lib('lib');
use ToyBox::LinearSMO;

my $smo = ToyBox::LinearSMO->new();

$smo->add_instance(attributes => {x => 1},  label => 1);
$smo->add_instance(attributes => {x => -1}, label => -1);
$smo->add_instance(attributes => {y => 1},  label => -1);
$smo->add_instance(attributes => {y => -1}, label => 1);

$smo->add_instance(attributes => {x => 2},  label => 1);
$smo->add_instance(attributes => {y => 2},  label => -1);

$smo->train(progress_cb => 'verbose');

print Dumper($smo);

print $smo->predict(attributes => {x => -0.5}), "\n";
print $smo->predict(attributes => {x => -0.5, y => -0.5}), "\n";

