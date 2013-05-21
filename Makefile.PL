use 5.006;
use strict;
use warnings FATAL => 'all';
use ExtUtils::MakeMaker;

sub MY::constants {
   package MY;   # To help SUPER work right 
   my $self = shift; 
   my @m; 

   push(@m, $self->SUPER::constants(@_)); 
   push(@m, 
	"PRGS =  baseml codeml basemlg mcmctree pamp evolver yn00 chi2\n",
	"LIBS = -lm # -lM\n",
	"PAMLDIR = paml4.7/\n",
	"SRC = \$(PAMLDIR)src/\n",
	"BASEML = \$(SRC)baseml.c\n",
	"CODEML = \$(SRC)codeml.c\n",
	"BASEMLG = \$(SRC)basemlg.c\n",
	"MCMCTREE = \$(SRC)mcmctree.c\n",
	"PAMP = \$(SRC)pamp.c\n",
	"EVOLVER = \$(SRC)evolver.c\n",
	"YN00 = \$(SRC)yn00.c\n",
	"CHI2 = \$(SRC)chi2.c\n",
	"TOOLS = \$(SRC)tools.c\n",
	"HEAD = \$(SRC)paml.h\n",
	"TREESUB = \$(SRC)treesub.c\n",
	"TREESPACE = \$(SRC)treespace.c\n");

   join "", @m; 
} 

WriteMakefile(
    NAME             => 'Mecom',
    AUTHOR           => q{Hector Valverde <hvalverde@uma.es>},
    VERSION_FROM     => 'lib/Mecom.pm',
    ABSTRACT_FROM    => 'lib/Mecom.pm',
    LICENSE          => 'Perl_5',
    EXE_FILES 	     => ['bin/mecom'],
    # Params needed to install PAML
    MAN3PODS         => {
	"bin/mecom" => "blib/man3/mecom.3pm"
    },
    MIN_PERL_VERSION => 5.006,
    CONFIGURE_REQUIRES => {
        'ExtUtils::MakeMaker' => 0,
    },
    PREREQ_PM => {
	"Bio::Perl" => "0",
	"Bio::Tools::Run::Phylo::PAML::Yn00" => "0",
	"Carp" => "0",
	"Getopt::Long" => "0",
	"Number::FormatEng" => "0",
	"Statistics::Basic" => "0",
	"Statistics::Zed" => "0",
	"IPC::Run" => "0",
    },
    BUILD_REQUIRES => {
        'Test::More' => 0,
    },
    dist  => { COMPRESS => 'gzip -9f', SUFFIX => 'gz', },
    clean => { FILES => 'Mecom-*'},
    realclean => { FILES => '$(PAMLDIR)bin/baseml
			    $(PAMLDIR)bin/codeml
			    $(PAMLDIR)bin/basemlg
			    $(PAMLDIR)bin/mcmctree
			    $(PAMLDIR)bin/pamp
			    $(PAMLDIR)bin/evolver
			    $(PAMLDIR)bin/chi2
			    $(PAMLDIR)bin/yn00'},
    depend => { '$(INST_SCRIPT)/mecom' => 'paml mvbin',
	        #pure_install => 'installpaml'
		},
);

package MY;

sub MY::postamble {'

# lets install paml4.7

paml : $(PRGS)

baseml : $(BASEML) $(TOOLS) $(TREESUB) $(TREESPACE) $(HEAD)
	$(CC) $(CFLAGS) -o $@ $(BASEML) $(TOOLS) $(LIBS)
basemlg : $(BASEMLG)  $(TOOLS) $(TREESUB) $(TREESPACE) $(HEAD) 
	$(CC) $(CFLAGS) -o $@ $(BASEMLG) $(TOOLS) $(LIBS)
codeml : $(CODEML)  $(TOOLS) $(TREESUB) $(TREESPACE) $(HEAD) 
	$(CC) $(CFLAGS) -o $@ $(CODEML) $(TOOLS) $(LIBS)
evolver : $(EVOLVER)  $(TOOLS) $(TREESUB) $(TREESPACE) $(HEAD) 
	$(CC) $(CFLAGS) -o $@ $(EVOLVER) $(TOOLS) $(LIBS)
pamp : $(PAMP)  $(TOOLS) $(TREESUB) $(TREESPACE) $(HEAD) 
	$(CC) $(CFLAGS) -o $@ $(PAMP) $(TOOLS) $(LIBS)
mcmctree : $(MCMCTREE)  $(TOOLS) $(TREESUB) $(TREESPACE) $(HEAD) 
	$(CC) $(CFLAGS) -o $@ $(MCMCTREE) $(TOOLS) $(LIBS)
yn00 : $(YN00) $(TOOLS) $(HEAD)
	$(CC) $(CFLAGS) -o $@ $(YN00) $(TOOLS) $(LIBS)
chi2 : $(CHI2)
	$(CC) $(CFLAGS) -o $@ $(CHI2) $(LIBS)
	
mvbin : $(PRGS)
	mv $(PRGS) $(PAMLDIR)bin/.
	
installpaml : all
	cp -R $(PAMLDIR) $(INSTALLBIN)/$(PAMLDIR)

 
'}