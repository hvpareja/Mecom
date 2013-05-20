use 5.006;
use strict;
use warnings FATAL => 'all';
use ExtUtils::MakeMaker;

WriteMakefile(
    NAME             => 'Mecom',
    AUTHOR           => q{Hector Valverde <hvalverde@uma.es>},
    VERSION_FROM     => 'lib/Mecom.pm',
    ABSTRACT_FROM    => 'lib/Mecom.pm',
    LICENSE          => 'Perl_5',
    PL_FILES         => {},
    EXE_FILES 	     => ['bin/mecom'],
    MAN3PODS         => {
	"bin/mecom" => "blib/man3/mecom.3pm"
    },
    MIN_PERL_VERSION => 5.006,
    CONFIGURE_REQUIRES => {
        'ExtUtils::MakeMaker' => 0,
    },
    PREREQ_PM => {
	"Bio::Structure::Model" => "0",
	"Bio::Structure::IO::pdb" => "0",
	"Bio::Structure::SetcStr::DSSP::Res" => "0",
	"Bio::SimpleAlign" => "0",
	"Bio::Align::Utilities" => "0",
	"Bio::Align::IO" => "0",
	"Bio::SeqIO" => "0",
	"Bio::Tools::Run::Phylo::PAML::Yn00" => "0",
	"Carp" => "0",
	"Getopt::Long" => "0",
	"Number::Formating" => "0",
	"Statistics::Basic" => "0",
	"Statistics::Zed" => "0",
	"IPC::Run" => "0",
    },
    BUILD_REQUIRES => {
        'Test::More' => 0,
    },
    dist  => { COMPRESS => 'gzip -9f', SUFFIX => 'gz', },
    clean => { FILES => 'Mecom-*'},
);
