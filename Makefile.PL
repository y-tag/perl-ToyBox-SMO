use inc::Module::Install;

name          'ToyBox-SMO';
all_from      'lib/ToyBox/SMO.pm';
license       'MIT';

test_requires 'Test::More';

WriteAll;
