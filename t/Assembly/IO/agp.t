use strict;
use warnings;
my %ASSEMBLY_TESTS;

BEGIN {
  use lib '.';
  use Bio::Root::Test;
  
  test_begin( -tests => 459+5,
	      # Not sure what to set here
	      -requires_modules => [ ],
	    );
  
  use_ok('Bio::Seq');
  use_ok('Bio::LocatableSeq');
  use_ok('Bio::Assembly::IO');
}

use Bio::Root::IO;

#my ($aio, $assembly, @contig_seq_ids, @singlet_ids, @contig_ids, @all_seq_ids);

my $file = 'chr_from_scaffold_WGS.agp';



# Test basic object construction
ok my $aio = Bio::Assembly::IO->
  new( -file   => test_input_file($file),
       -format => 'agp' ),
  "init agp IO object";
isa_ok($aio, 'Bio::Assembly::IO');



# Test getting assemblies out of the object
while (my $assembly = $aio->next_assembly) {
  isa_ok($assembly, 'Bio::Assembly::Scaffold');
}



# Test reopening the file
ok $aio = Bio::Assembly::IO->
  new( -file => test_input_file($file),
       -format => 'agp' ), "reopen";
isa_ok($aio, 'Bio::Assembly::IO');



# Look at a returned assembly a bit
my $assembly = $aio->next_assembly;

isa_ok($assembly, 'Bio::Assembly::Scaffold');

is ($assembly->id, 'chrY');

is ($assembly->get_nof_contigs, 23);
is ($assembly->get_nof_contig_seqs, 23);
is ($assembly->get_nof_singlets, 0);
is ($assembly->get_all_seq_ids, 23);
is ($assembly->get_nof_seqs, 23);
is ($assembly->get_contig_seq_ids, 23);
is ($assembly->get_contig_ids, 23);

is ($assembly->get_singlet_ids, 0);





# Look at assembly annotations
my $assembly_annot = $assembly->annotation;
isa_ok($assembly_annot, 'Bio::Annotation::Collection');

foreach my $key ( $assembly_annot->get_all_annotation_keys() ) {
  my @values = $assembly_annot->get_Annotations($key);
  foreach my $value ( @values ) {
    # value is an Bio::AnnotationI, and defines a "as_text" method
    print "Annotation ", $key, " stringified value ", $value->as_text, "\n";
  }
}

__END__
foreach my contig ( $assembly->_contig) {
  isa_ok($contig, 'Bio::Assembly::Contig');
}


__END__


ok $assembly = $aio->next_assembly, "get agp assy";
is( $assembly->get_nof_contigs, 23, "got all contigs");


ok(@contig_seq_ids = $assembly->get_contig_seq_ids, "get_contig_seq_ids");
is(@contig_seq_ids, 334);
