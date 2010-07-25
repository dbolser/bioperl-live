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



# Look at a contig

is ($assembly->get_contig_by_id( 'foo' ), undef );

my $contig = $assembly->get_contig_by_id( 'EG1_scaffold1' );

isa_ok($contig, 'Bio::Assembly::Contig');

$contig = ($assembly->remove_contigs_by_id( 'EG1_scaffold1' ))[0];

isa_ok($contig, 'Bio::Assembly::Contig');

is ($assembly->get_nof_contigs, 22);
is ($assembly->get_nof_contig_seqs, 22);
is ($assembly->get_nof_singlets, 0);
is ($assembly->get_all_seq_ids, 22);
is ($assembly->get_nof_seqs, 22);
is ($assembly->get_contig_seq_ids, 22);
is ($assembly->get_contig_ids, 22);
is ($assembly->get_singlet_ids, 0);


## Look at the contig

is ($contig->source, 'not specified');

isa_ok($contig->assembly, 'Bio::Assembly::Scaffold');

is ($contig->assembly->id, 'chrY');

is($contig->strand, 0);

print $contig->upstream_neighbor, "\n";
print $contig->downstream_neighbor, "\n";


# Look at contig features
my $contig_features = $contig->get_features_collection;
isa_ok($contig_features, 'Bio::SeqFeature::Collection');

is ($contig_features->feature_count, 2);

foreach my $feat ( $contig_features->get_all_features ) {
  isa_ok($feat, 'Bio::SeqFeatureI');
  
  print
    "Feature from ", $feat->start, " to ", $feat->end,
      " Primary tag  ", $feat->primary_tag,
	", produced by ", $feat->source_tag(), "\n";
}






__END__



is ($assembly->all_contigs, 22);

foreach my $contig ( $assembly->all_contigs) {
  isa_ok($contig, 'Bio::Assembly::Contig');
}







ok $assembly = $aio->next_assembly, "get agp assy";
is( $assembly->get_nof_contigs, 23, "got all contigs");


