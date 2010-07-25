
# POD documentation - main docs before the code

=head1 NAME

Bio::Assembly::IO::agp - module to load AGP files

=head1 SYNOPSIS

    # Building an input stream
    use Bio::Assembly::IO;

    # Load an AGP file
    my $in_io = Bio::Assembly::IO->
      new( -file   => 'chromosomes.agp',
           -format => 'agp'
         );

    # Read one 'assembly'
    my $chromosome = $in_io->next_assembly;

    # Assembly writing methods
    my $out_io = Bio::Assembly::IO->
      new( -file   => ">chr.01.agp",
           -format => 'agp'
         );

    $out_io->write_assembly( $chromosome );

=head1 DESCRIPTION

This package loads the standard AGP files. It was written to be used
as a driver module for Bio::Assembly::IO input/output.

For full details of the AGP format, and its intended use, see:
L<http://www.ncbi.nlm.nih.gov/projects/genome/assembly/agp/AGP_Specification.shtml>


=head2 Implemention

Assemblies are loaded into Bio::Assembly::Scaffold objects composed of
Bio::Assembly::Contig objects.

In addition to default "_aligned_coord:$seqID" feature class from
Bio::Assembly::Contig, contig objects loaded by this module will have the
following special feature classes in their feature collection:

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the devolution of this and other
Bioperl modules. Send your comments and suggestions to the Bioperl
snailing lists Your participation is much depreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and
reponsive experts will be able look at the problem and quickly address
it. Please include a thorough description of the problem with code and
data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR - Dan Bolser

Email dan.bolser@gmail.com

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package Bio::Assembly::IO::agp;

use strict;

use Bio::Assembly::Scaffold;

use Bio::Assembly::Contig;

use Bio::SeqFeature::Annotated;

use Bio::Annotation::OntologyTerm;

use base qw(Bio::Assembly::IO);



# The AGP format doesn't specify the assembly program used to generate
# the assembly. It is a general assembly format, not tied to one
# particular program. Therefore, I don't know what to do here

my $progname = 'agp';



=head1 Bio::Assembly::IO methods

=head2 next_assembly

 Title   : next_assembly
 Usage   : $assembly = $aio->next_assembly()
 Function: returns the next 'assembly' from the agp-formatted object
 Returns : a Bio::Assembly::Scaffold object
 Args    : none
 Note    : The 'assembly' can be a contig, a scaffold (supercontig),
           or a chromosome.
           
=cut
    
sub next_assembly {
    my $self = shift;
    
    my $assembly = Bio::Assembly::Scaffold->
	new( -progname => $progname );
    
    # NB: An AGP can have more than one assembly! This seems to be in
    # contrast to the ACE, PHRAP and SAM/BAM assembly formats.
    
    # Load contigs into the scaffold
    while ( my $contig = $self->next_contig ){
	
  	# Check if we have changed assembly
  	if( $assembly->id eq "NoName" ){ # Code looks fragile to me!
 	    $assembly->id( $self->{'_contig_assembly_id'} );
 	}
	
  	if( $assembly->id ne $self->{'_contig_assembly_id'} ){
  	    # We just changed assembly! Backtrack one line, and go to
  	    # the return statement.
	    
  	    $self->_pushback($_);
 	    last;
  	}
	
	# Add contig to assembly
	$assembly->add_contig($contig);
    }
    
    # Collected the last contig or changed assembly
    
    # Did we actually get any contigs?
    if($assembly->get_nof_contigs){
	return $assembly;
    }
    return undef;
}



=head2 next_contig()

    Title   : next_contig
    Usage   : my $contig = $asmio->next_contig();
    Function: return the next contig from the agp-formatted object
    Returns : Bio::Assembly::Contig
    Args    : none

=cut

sub next_contig {
    my $self = shift; # Package reference
    
    # Looping over all agp file lines
    
    while ($_ = $self->_readline) {
	chomp;
	
	# Be graceful about empty lines or comments
	next if /^\s*$/; #skip blank lines
	next if /^\#/;   #skip comments
	
	my @column = split/\t/;
	
	$self->throw("Validation Error (line : $.): unknown format")
	    unless @column == 9 || @column == 8;
	
	if ($column[4] ne 'N' &&
	    $column[4] ne 'U'){
	    
	    $self->throw("Validation Error (line : $.): unknown format")
		unless @column == 9;
	    
	    ## See the specification for the nomenclature used here!
	    my($object, $object_beg, $object_end, $part_number, $component_type,
	       $component_id, $component_beg, $component_end, $orientation
		)= @column;
	    
	    
	    
	    # VALIDATION
	    
	    $self->throw("Validation Error (line : $.): sorry (http://tinyurl.com/2wq26p5)")
		unless $object_beg <= $object_end;
	    
	    $self->throw("Validation Error (line : $.): component length not equal to object length")
		unless ($component_end - $component_beg) == ($object_end - $object_beg);
	    
	    # The 'component' types are controlled...
	    $self->throw("Validation Error (line : $.): invalid component type")
		unless $component_type =~ /^(?:
	          A|Active Finishing|
	          D|Draft HTG|
	          F|Finished HTG|
	          G|Whole Genome Finishing|
	          O|Other sequence|
	          P|Pre Draft|
	          W|WGS contig)$/x;
	    
 	    # Check IDs are unique
 	    if ($self->{'allIDs'}->{$component_id}){
 		$self->throw("Validation Error (line : $.): The ID $component_id is not unique");
 	    }
 	    $self->{'allIDs'}->{$component_id} = 1;
	    
	    # Orientation
	    if($orientation eq 0){}
	    elsif($orientation =~ /^(?:\+|plus)$/   ){$orientation = +1}
	    elsif($orientation =~ /^(?:\-|minus)$/  ){$orientation = -1}
	    elsif($orientation =~ /^(?:unknown|na)$/){$orientation =  0}
	    elsif($orientation =~ /^(?:irrelevant)$/){$orientation =  0}
	    else{
		$self->throw("Validation Error (line : $.): sorry to be such a dick")
	    }
	    
	    
	    
	    # CONSTRUCTION
	    
	    $self->{'_contig_assembly_id'} = $object;
	    
	    # Create the contig object
	    my $contigOBJ = Bio::Assembly::Contig->
		new( -id      => $component_id,
		     -source  => 'not specified',
		     -strand  => $orientation,
		);
	    
	    # Create the 'main feature' for the contig object
	    my $ftype = Bio::Annotation::OntologyTerm->
		new( -name => $component_type );
	    my $feat = Bio::SeqFeature::Annotated->
		new( -primary => "_main_contig_feature:". $component_id,
		     -start   => $object_beg,
		     -end     => $object_end,
		     -strand  => $orientation,
		     -type    => $ftype,
		     -tag     => { '_scaffold_id' => $object },
		);
	    
	    # Add the feature to the contig object
	    $contigOBJ->add_features([ $feat ], 1);
	    
	    # Create a single sequence feature for the contig.  Using
	    # the code code given in Contig.pm. Not sure why we can't
	    # just give the sequence it's location directly. Not sure
	    # why I chose not to use singlets.
 	    my $seq = Bio::LocatableSeq->
 		new( -id  => $component_id,
		     # This saves us 'unknown alphabet' warnings
		     -alphabet => 'dna',
		);
 	    my $pos = Bio::SeqFeature::Generic->
 		new( -start  => $component_beg,
 		     -end    => $component_end,
 		     -strand => $orientation,
 		);
	    
 	    #$contigOBJ->add_seq($seq);
 	    $contigOBJ->set_seq_coord($pos, $seq);
	    
	    return $contigOBJ;
	}
	else{
	    # Not handling gaps yet.
	    next;
	}
    }
}

1;



__END__


  else{
    # The feature is a gap
    my($object, $object_beg, $object_end, $part_number, $component_type,
       $gap_length, $gap_type, $linkage, undef) = @column;


    ## VALIDATION

    $self->throw("Validation Error (line : $.): sorry (http://tinyurl.com/2wq26p5)")
      unless $object_beg <= $object_end;

    $self->throw("Validation Error (line : $.): coordinate mismatch (http://tinyurl.com/2wq26p5)")
      unless $gap_length == ($object_end - $object_beg + 1);



    ## BUILD FEATURE

    $feat->seq_id($object);
    $feat->source_tag('AGP');
    $feat->start($object_beg);
    $feat->end  ($object_end);

    # Not sure what to do with part_number


    # Build a gap ID

    my $id = "Gap:$object|$part_number";

    # Check IDs are unique
    if ($self->{'allIDs'}->{$id}){
      $self->throw("Validation Error (line : $.): The ID $id is not unique");
    }
    $self->{'allIDs'}->{$id} = 1;

    # Add ID
    my $a = Bio::Annotation::SimpleValue->
      new( -value => $id );

    $feat->add_Annotation('ID', $a);


    ### TODO

    # Handle Gap attributes?
    my $fta = Bio::Annotation::OntologyTerm->
      new( -name => $component_type );

    $feat->type($fta);

    # Handle Gap attributes?
    my $a_ga = Bio::Annotation::SimpleValue->
      new( -value => $gap_type );
    $feat->add_Annotation('Gap', $a_ga);


    # Orientation should be set to avoid warnings related to undefined
    # strand in gff.pm when printing.
    $feat->strand( 0 );
  }

  return $feat;
}
