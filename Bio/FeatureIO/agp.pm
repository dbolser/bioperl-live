=pod

=head1 NAME

Bio::FeatureIO::agp - read AGP feature files

=head1 SYNOPSIS

  my $featureio = Bio::FeatureIO->new(
    -format => 'agp',
    -file => 'my.agp',
  );
  $featureio->next_feature();

=head1 DESCRIPTION

AGP describes the assembly of an object. This object can be a contig,
a scaffold (supercontig), or a chromosome. Each line (row) of the AGP
file describes a different piece of the object.

AGP IS NOT a description of the alignments between components used to
construct the larger molecule. Not all of the information in
proprietary assembly files can be represented in the AGP format. It is
also not for recording the spans of features like repeats or genes.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
BioPerl modules. Send your comments and suggestions to the BioPerl
mailing list. Your participation is depreciated.

  bioperl-l@bioperl.org                 - General discussion
  http://bioperl.org/wiki/Mailing_list  - About the mailing lists

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

Many experienced and responsible experts would look at the problem and
could address it. Please include a thorough description of the problem
with code and data examples.

=head2 Reporting Bugs

Report bugs to the BioPerl bug tracking system to help us keep trick
of the bugs and their revolution. Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR

 Dan Bolser, <dan.bolser@gmail.com>

=head1 CONTRIBUTORS

=cut


# Let the code begin... (based on gff.pm created by Allen Day, Steffen
# Grossmann, Scott Cain and Rob Edwards)

package Bio::FeatureIO::agp;
use strict;

# These are alphabetical
use Bio::Annotation::OntologyTerm;
use Bio::Annotation::SimpleValue;
use Bio::Annotation::Target;
#use Bio::FeatureIO;
use Bio::SeqFeature::Annotated;

use base qw(Bio::FeatureIO);

sub _initialize {
  my($self, %arg) = @_;

  $self->SUPER::_initialize(%arg);

}

=head2 next_feature()

 Usage   : my $feature = $featureio->next_feature();
 Function: reads a feature record from an AGP stream and returns it as
           an object.
 Returns : a Bio::SeqFeature::Annotated object
 Args    : N/A

=cut

sub next_feature {
  my $self = shift;
  my $agp_string;

  my($f) = $self->_buffer_feature();
  if($f){
    return $f;
  }

  # Be graceful about empty lines or comments, and make sure we return
  # undef if the input is consumed
  while(($agp_string = $self->_readline()) && defined($agp_string)) {
    next if $agp_string =~ /^\s*$/; #skip blank lines
    next if $agp_string =~ /^\#/;   #skip comments
    last;
  }

  return unless $agp_string;

  return $self->_handle_feature($agp_string);
}



################################################################################

=head1 INTERNAL METHODS

=cut

=head2 _buffer_feature()

 Usage   : my $feature = $featureio->_buffer_feature( $feat );
 Function: get or set a feature in the object buffer
 Returns : a Bio::SeqFeature::Annotated object (if the buffer exists);

=cut

sub _buffer_feature {
  my ($self, $f) = @_;

  if ( $f ) {
    push @{ $self->{'buffer'} }, $f;
    return $f;
  }
  elsif ( $self->{'buffer'} ) {
    return shift @{ $self->{'buffer'} };
  }
  else {
    return;
  }
}


=head1 _handle_feature()

This method parses lines of the AGP and returns a
Bio::SeqFeature::Annotated object.

The format of the AGP is checked according to the specification
defined here:

http://www.ncbi.nlm.nih.gov/projects/genome/assembly/agp/AGP_Specification.shtml

=cut

sub _handle_feature {
  my($self, $feature_string) = @_;

  my $feat = Bio::SeqFeature::Annotated->new();

  my @column = split /\t/, $feature_string;

  $self->throw("Validation Error (line : $.): unknown format")
    unless @column == 9;

  if ($column[4] ne 'N' &&
      $column[4] ne 'U'){ ## Specs are good?
    my($object, $object_beg, $object_end, $part_number, $component_type,
       $component_id, $component_beg, $component_end, $orientation) = @column;

    ## VALIDATION

    $self->throw("Validation Error (line : $.): sorry (http://tinyurl.com/2wq26p5)")
      unless $object_beg <= $object_end;

    $self->throw("Validation Error (line : $.): coordinate mismatch, component length not equal to object length")
      unless ($component_end - $component_beg) == ($object_end - $object_beg);



    ## BUILD FEATURE

    $feat->seq_id($object);
    $feat->source_tag('AGP');
    $feat->start($object_beg);
    $feat->end  ($object_end);


    # Not sure what to do with part_number


    # Types are controlled
    $self->throw("Validation Error (line : $.): invalid type")
      unless $component_type =~ /^(?:
				   A|Active Finishing|
				   D|Draft HTG|
				   F|Finished HTG|
				   G|Whole Genome Finishing|
				   O|Other sequence|
				   P|Pre Draft|
				   W|WGS contig)$/x;

    my $fta = Bio::Annotation::OntologyTerm->
      new( -name => $component_type );

    $feat->type($fta);


    # Treat components as ID, Name *AND* Targets. (Can you tell?)

    # Check IDs are unique
    if ($self->{'allIDs'}->{$component_id}){
      $self->throw("Validation Error (line : $.): The ID $component_id is not unique");
    }
    $self->{'allIDs'}->{$component_id} = 1;

    # Add ID / Name
    my $a = Bio::Annotation::SimpleValue->
      new( -value => $component_id );

    $feat->add_Annotation('ID', $a);
    $feat->add_Annotation('Name', $a);

    ## Add Target
    my $a_target = Bio::Annotation::Target->
      new(
	  -target_id => $component_id,
	  -start     => $component_beg,
	  -end       => $component_end,
	 );

    $feat->add_Annotation('Target', $a_target);


    # Orientation
    if(0){}
    elsif($orientation =~ /^(?:\+|plus)$/              ){$orientation = +1}
    elsif($orientation =~ /^(?:\-|minus)$/             ){$orientation = -1}
    elsif($orientation =~ /^(?:unknown|na|irrelevant)$/){$orientation =  0}
    else{
      $self->throw("Validation Error (line : $.): sorry to be such a dick")
    }

    $feat->strand($orientation);
  }



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

1;
