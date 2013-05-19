package Mecom;
# -----------------------------------------------------------------------------
# Molecular Evolution of Protein Complexes Contact Interfaces
# -----------------------------------------------------------------------------
# @Authors:  HŽctor Valverde <hvalverde@uma.es> and Juan Carlos Aledo
# @Date:     May-2013
# @Location: Depto. Biolog’a Molecular y Bioqu’mica
#            Facultad de Ciencias. Universidad de M‡laga
#
# Copyright 2013 Hector Valverde and Juan Carlos Aledo.
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of either: the GNU General Public License as published
# by the Free Software Foundation; or the Artistic License.
#
# See http://dev.perl.org/licenses/ for more information.
# -----------------------------------------------------------------------------

use 5.006;
use strict;
no strict "refs";
our $VERSION = '1.06';

# Own modules
use Mecom::Contact; 
use Mecom::Surface;
use Mecom::Subsets;
use Mecom::Report;
use Mecom::Align::Subset;         # Avaliable from CPAN
use Mecom::EasyYang;              # Avaliable from CPAN
use Mecom::Statistics::RatioVariance;  # Avaliable from CPAN

# External modules
use Bio::Structure::Model;
use Bio::Structure::IO::pdb;
use Bio::SimpleAlign;
use Bio::Align::Utilities qw(:all);
use warnings;
use Carp;
use Bio::AlignIO;

# -----------------------------------------------------------------------------
# Class data                                                          Chap.  1
# -----------------------------------------------------------------------------
{
    # A list of all attributes wiht default values and read/write/required
    # properties
    my %_attribute_properties = (
        
        # Required (pdb become optional if contactfile is set)
        _pdb                   => ["-"        , ""         ],
        _contactfile           => ["-"        , ""         ],
        _alignment             => [""         , "required" ],
        _chain                 => [""         , "required" ],
        
        # Optional
        _pth                   => [4          , ""         ],
        _sth                   => [0.05       , ""         ],
        _sthmargin             => [0          , ""         ],
        _contactwith           => ["-"        , ""         ],
        _informat              => ["fasta"    , ""         ],
        _oformat               => ["clustalw" , ""         ],
        _gc                    => [0          , ""         ],
        _ocontact              => ["data.str" , ""         ],
        _dsspbin               => ["dssp"     , ""         ],
        _report                => ["r.html"   , ""         ],
        
        # Processed information
        _struct_data           => [[]         , ""         ],
        _lists                 => [{}         , ""         ],
        _sub_alns              => [{}         , ""         ],
        _paml_res              => [{}         , ""         ],
        _stats                 => [{}         , ""         ],
        
    );
    
    # Global variable to keep count of existing objects
    my $_count = 0;
    # The list of all attributes
    sub _all_attributes {
        keys %_attribute_properties;
    }
    # Check if a given property is set for a given attribute
    sub _permissions{
        my ($self,$attribute, $permissions) = @_;
        $_attribute_properties{$attribute}[1] =~ /$permissions/;
    }
    # Return the default value for a given attribute
    sub _attribute_default{
        my ($self,$attribute) = @_;
        $_attribute_properties{$attribute}[0];
    }
    # Manage the count of existing objects
    sub get_count{
        $_count;
    }
    sub _incr_count{
        ++$_count;
    }
    sub _decr_count{
        --$_count;
    }

}

    

# -----------------------------------------------------------------------------
#                                                                            1.

#                                   #---#                                    #


# -----------------------------------------------------------------------------
# Constructor                                                         Chap.  2
# -----------------------------------------------------------------------------
sub new{
    
        
    my ($class, %arg) = @_;
    my $self = bless {}, $class;

    foreach my $attribute ($self->_all_attributes()){
        
        # E.g. attribute = "_name", argument = "name"
        my ($argument) = ($attribute =~ /^_(.*)/);
        
        # If explicitly given
        if($arg{$argument}){
            $self->{$attribute} = $arg{$argument};
        }
        
        # If not given but required
        elsif($self->_permissions($attribute,'required')){
            croak("No $argument specified as required");
        }
        
        # Set to default
        else{
            $self->{$attribute} = $self->_attribute_default($attribute);            
        }
        
    }
    
    # Test if DSSP and PAML exists
    if( -e $self->get_dsspbin ){
        croak("The DSSP program (".$self->get_dsspbin.") is missing.
               See documentation to solve this error.");
    }
    if($ENV{PAMLDIR} eq ""){
        croak("PAML software must be correctly installed in your system.
               See documentation to solve this error.");
    }
    
    # Called $class because it is a gobal method
    $class->_incr_count;
    

    # Input requirement
    if($self->get_pdbfilepath eq "-" && $self->get_contactfilepath eq "-"){
        croak("Input file (pdb or contact) is required")
    }
    
    return $self;
    
}
# -----------------------------------------------------------------------------
#                                                                            2.

#                                   #---#                                    #

# -----------------------------------------------------------------------------
# Methods                                                             Chap.  3
# -----------------------------------------------------------------------------
# Run complete
sub run{
    
    my $self = $_[0];
    # Structural data
    $self->run_struct;
    
    # Filtering
    $self->run_filtering;
    
    # Sub-alignments
    $self->run_subalign;
    
    # Yang
    $self->run_yang;
    
    # Stats
    $self->run_stats1;
    
}
# Get Structural data
sub run_struct{
    
    my $self = $_[0];
    if($self->get_contactfilepath eq "-"){
    
        # Contacts
        my @contact_data = $self->_run_contacts;
        # Surface
        $self->set_structdata($self->_run_surface(@contact_data));
    
    }else{
    
        $self->set_structdata($self->_run_contactFromFile);
        $self->set_ocontact($self->get_contactfilepath);
    
    }
    
    return 1;
    
}
# Get contacts
sub _run_contacts{
    
    my ($self) = $_[0];
    my @contact_data;
    my $stream = Bio::Structure::IO->new(-file => $self->get_pdbfilepath,
                                         -format => 'PDB')
                                 or die "\nInvalid pdb file.\n";

    # Switch if contact file or pdbfile
    if($self->get_contactfilepath() eq "-"){
   
    print "\tCalculating contacts ... it may take a few minutes, please wait.\n";
                                 
    my $obj_contact = new Mecom::Contact("th"         => $self->get_pth,
                                               "pdb"        => $stream,
                                               "chain"      => $self->get_chain);
    
    my $residue_contacts = $obj_contact->contacts;
    @contact_data = @$residue_contacts;
    
    return @contact_data;
    
    }
}
# Get surface
sub _run_surface{
    
    my ($self,@contact_data) = @_;
    my @structural_info;
    
    my $plus_dssp  = Mecom::Surface::dssp($self->get_pdbfilepath,
                                                $self->get_chain,
                                                $self->get_sth,
                                                $self->get_sthmargin,
                                                $self->get_dsspbin,
                                                @contact_data);
    # A new asignation just for clarity
    my @dssp = @$plus_dssp;
    @structural_info = @dssp;
    
    # Save this data (save time for further analisys)
    open OCONTACT, ">".$self->get_ocontact();
    foreach my $line (@dssp){
        print OCONTACT $line."\n";
    }
    close OCONTACT;
    
    return @structural_info
    
}
# Get structural info from contact file
sub _run_contactFromFile{
    
    my $self = $_[0];
    my $contact_file = $self->get_contactfilepath;
    
    print "\tGetting data from contact file\n";
    open CONTACT, $contact_file;
    my @structural_info = <CONTACT>;
    close CONTACT;
    
    return @structural_info;
    
}
# Get subsets
sub run_filtering{
    
    my $self = $_[0];
    
    # With the help of regular expresions, the program will extract a residue list
    # to build a new alignment
    my %subsets_list = Mecom::Subsets->build($self->get_chain,
                                               $self->get_contactwith,
                                               @{$self->get_structdata});
    
    # If a subset is empty, it will be deleted
    foreach my $key (keys %subsets_list){
        if(!$subsets_list{$key}[0]){ delete $subsets_list{$key}; }
    }
    
    $self->set_lists(%subsets_list);
    return 1;
    
}
# Get sub-alignments
sub run_subalign{
    
    my $self = $_[0];
    my %lists = %{$self->get_lists};
    my $obj = Mecom::Align::Subset->new(file   => $self->get_alignfilepath,
                                      format => $self->get_informat);
    
    my %aln_subsets;
    foreach my $key (keys %lists){
        # This function returns a Bio::SimpleAlign object
        $aln_subsets{$key} = $obj->build_subset($lists{$key});
    }

    $self->set_subalns(%aln_subsets);
    return 1;
    
}
# Get evolutionary analisys
sub run_yang{
    
    my $self = $_[0];
    my %alns = %{$self->get_subalns};
    my %paml_results;
    foreach my $aln_key (keys %alns){
        $paml_results{$aln_key} = {Mecom::EasyYang->
                                   yang($alns{$aln_key},$self->get_gc)};
    }
    
    $self->set_pamlres(%paml_results);
    return 1;
    
}
# Get statistics 1
sub run_stats1{
    
    my $self = $_[0];
    my %paml_results = %{$self->get_pamlres};
    my %results;
    foreach my $category (keys %paml_results){
        foreach my $other (keys %paml_results){
            if($other ne $category){
                
                my ($x,$y,$var_x,$var_y) = ($paml_results{$category}{dN},
                                            $paml_results{$other}{dN},
                                            $paml_results{$category}{dN_VAR},
                                            $paml_results{$other}{dN_VAR});
                my @x = @$x;
                my @y = @$y;
                my @var_x = @$var_x;
                my @var_y = @$var_y;
                if($#x == $#y){
                #print $#x." -- ".$#y."\n";
                    $results{$category." vs ".$other} =
                                        {Mecom::Statistics::RatioVariance->calc($x,
                                                                        $y,
                                                                        $var_x,
                                                                        $var_y
                                                                        )};
                }
            }
        }
    }
    
    # Stats report
    foreach my $key (keys %results){
        if(!$results{$key}{standar_deviation}){ # Sets are equal
            delete $results{$key};
        }
        
    }
    
    $self->set_stats(%results);
    return 1;
    
}
# Get report
sub run_report{
    
    my $self = $_[0];
    my $html_report;
    my @input_info = (
                  $self->get_pdbfilepath,
                  $self->get_contactfilepath,
                  $self->get_alignfilepath,
                  $self->get_chain,
                  $self->get_gc,
                  $self->get_contactwith,
                  $self->get_pth,
                  $self->get_sth,
                  $self->get_sthmargin,
                  $self->get_informat,
                  $self->get_oformat,
                  $self->get_ocontact
                  );
    
    $html_report = Mecom::Report->input_information  (@input_info);
    $html_report.= Mecom::Report->struct_information ($self->get_pdbfilepath,
                                                            $self->get_structdata,
                                                            $self->get_chain,
                                                            $self->get_ocontact);
    $html_report.= Mecom::Report->codon_lists        ($self->get_lists);  
    $html_report.= Mecom::Report->sub_alignments     ($self->get_chain,
                                                        $self->get_oformat,
                                                        $self->get_subalns);  
    $html_report.= Mecom::Report->yang_report        ($self->get_pamlres);     
    $html_report.= Mecom::Report->stats1             ($self->get_stats);           
    
    open HTML, ">".$self->get_report;
    print HTML $html_report;
    close HTML,
    
    return 1;
    
}
# -----------------------------------------------------------------------------
#                                                                            3.

#                                   #---#                                    #

# -----------------------------------------------------------------------------
# Auxiliar Methods                                                    Chap.  4
# -----------------------------------------------------------------------------
# Concatenate alignments
# This is a global function and must be called as
# Mecom::Complex->cat_aln(@alns);
sub cat_aln{
    
    # Bio::SimpleAlign objects as arguments
    my ($self, @alns) = @_;
    # Call function Bio::Align::Utilities->cat();
    my $merge_aln = cat(@alns);
    
    return $merge_aln;
    
}
# -----------------------------------------------------------------------------
#                                                                            4.

#                                   #---#                                    #

# -----------------------------------------------------------------------------
# Accesor Methods                                                    Chap.  5
# -----------------------------------------------------------------------------
# This kind of method is called Accesor
# Method. It returns the value of a key
# and avoid the direct acces to the inner
# value of $obj->{_file}.
sub get_pdbfilepath        { $_[0] -> {_pdb}         }
sub get_contactfilepath    { $_[0] -> {_contactfile} }
sub get_alignfilepath      { $_[0] -> {_alignment}   }
sub get_chain              { $_[0] -> {_chain}       }
sub get_pth                { $_[0] -> {_pth}         }
sub get_sth                { $_[0] -> {_sth}         }
sub get_sthmargin          { $_[0] -> {_sthmargin}   }
sub get_contactwith        { $_[0] -> {_contactwith} }
sub get_informat           { $_[0] -> {_informat}    }
sub get_oformat            { $_[0] -> {_oformat}     }
sub get_gc                 { $_[0] -> {_gc}          }
sub get_ocontact           { $_[0] -> {_ocontact}    }
sub get_dsspbin            { $_[0] -> {_dsspbin}     }
sub get_report             { $_[0] -> {_report}      }
# Proc
sub get_structdata         { $_[0] -> {_struct_data} }
sub get_lists              { $_[0] -> {_lists}       }
sub get_subalns            { $_[0] -> {_sub_alns}    }
sub get_pamlres            { $_[0] -> {_paml_res}    }
sub get_stats              { $_[0] -> {_stats}       }
# -----------------------------------------------------------------------------
#                                                                            5.

#                                   #---#                                    #

# -----------------------------------------------------------------------------
# Mutator Methods                                                     Chap.  6
# -----------------------------------------------------------------------------
sub set_pdbfilepath        { my ($self, $var) = @_;
                            $self->{_pdb} = $var if $var; }
sub set_contactfilepath    { my ($self, $var) = @_;
                            $self->{_contactfile} = $var if $var; }
sub set_alignfilepath      { my ($self, $var) = @_;
                            $self->{_alignment} = $var if $var;   }
sub set_chain              { my ($self, $var) = @_;
                            $self->{_chain} = $var if $var;       }
sub set_pth                { my ($self, $var) = @_;
                            $self->{_pth} = $var if $var;         }
sub set_sth                { my ($self, $var) = @_;
                            $self->{_sth} = $var if $var;         }
sub set_sthmargin          { my ($self, $var) = @_;
                            $self->{_sthmargin} = $var if $var;   }
sub set_contactwith        { my ($self, $var) = @_;
                            $self->{_contactwith} = $var if $var; }
sub set_informat           { my ($self, $var) = @_;
                            $self->{_informat} = $var if $var;    }
sub set_oformat            { my ($self, $var) = @_;
                            $self->{_oformat} = $var if $var;     }
sub set_gc                 { my ($self, $var) = @_;
                            $self->{_gc} = $var if $var;          }
sub set_ocontact           { my ($self, $var) = @_;
                            $self->{_ocontact} = $var if $var;    }
sub set_dsspbin            { my ($self, $var) = @_;
                            $self->{_dsspbin} = $var if $var;     }
sub set_report             { my ($self, $var) = @_;
                            $self->{_report} = $var if $var;     }
# Proc
sub set_structdata         { my ($self, @var) = @_;
                            $self->{_struct_data} = \@var if @var;}
sub set_lists              { my ($self, %var) = @_;
                            $self->{_lists} = \%var if %var;     }
sub set_subalns            { my ($self, %var) = @_;
                            $self->{_sub_alns} = \%var if %var;     }
sub set_pamlres            { my ($self, %var) = @_;
                            $self->{_paml_res} = \%var if %var;     }
sub set_stats              { my ($self, %var) = @_;
                            $self->{_stats} = \%var if %var;     }
# -----------------------------------------------------------------------------
#                                                                            6.

#                                   #---#                                    #

# -----------------------------------------------------------------------------
#                                                                     Chap.  7
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
#                                                                            7.

#                                   #---#                                    #

1;