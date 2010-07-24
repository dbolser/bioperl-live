use strict;
use warnings;
my %ASSEMBLY_TESTS;

BEGIN {
    use lib '.';
    use Bio::Root::Test;

    test_begin( -tests => 459+5,
                -requires_modules => [
                   ],
               );

    use_ok('Bio::Seq');
    use_ok('Bio::LocatableSeq');
    use_ok('Bio::Seq::Quality');
    use_ok('Bio::Assembly::IO');
    use_ok('Bio::Assembly::Singlet');
}

use Bio::Root::IO;

my ($aio, $assembly, @contig_seq_ids, @singlet_ids, @contig_ids, @all_seq_ids);
my $file = 'chr_from_contig_BAC.agp';

ok $aio = Bio::Assembly::IO->
  new( -file => test_input_file($file),
       -format => 'agp' ), "init agp IO object";
isa_ok($aio, 'Bio::Assembly::IO');
