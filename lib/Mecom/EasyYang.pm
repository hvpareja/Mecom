package Mecom::EasyYang;

# Library to implement the Yang ML algorithm
# This is just a PAML Wrapper
use Bio::Tools::Run::Phylo::PAML::Yn00;


sub yang {

    my ($self,$aln,$gen_code) = @_;

    #
    # Running Yang (Yn00)
    #
    # Create the Yn00 object
    my $yang = Bio::Tools::Run::Phylo::PAML::Yn00->new();
    
    # Set the genetic code
    $yang->set_parameter("icode", $gen_code);   
    
    # Append the alignment to Yn00 object
        $yang->alignment($aln);
    
    # Run Yang
        my ($rc,$parser) = $yang->run();
    
      
    my @dN     = ();
    my @dS     = ();
    my @dN_VAR = ();
    my @dS_VAR = ();
    my @N      = ();
    my @S      = ();
    my @omega  = ();
    my @kappa  = ();
    my @t      = ();
    
    # Take the ML Matrix
    my $out_str = "";
    $out_str   .= "----------------------------------------------------------------------\n";
    $out_str   .= "N\tS\tdN\tdN_SE\tdS\tdS_SE\tOmega\tKappa\tt\n";
    $out_str   .= "----------------------------------------------------------------------\n";
        
        while( my $result = $parser->next_result){
            
            my $MLmatrix = $result->get_MLmatrix();
            
            for(my $i=0;$i<=$#{$MLmatrix};$i++){
                
                for my $item (@{$MLmatrix}[$i]){
                    
                    for(my $j=$i+1;$j<=$#{$item};$j++){
                        
                        my $N     = $item->[$j]->{N};
                        my $S     = $item->[$j]->{S};
                        my $dN    = $item->[$j]->{dN};
                        my $dS    = $item->[$j]->{dS};
                        my $dN_SE = $item->[$j]->{dN_SE};
                        my $dS_SE = $item->[$j]->{dS_SE};
                        my $omega = $item->[$j]->{omega};
                        my $kappa = $item->[$j]->{kappa};
                        my $t     = $item->[$j]->{t};
                        
                        $dN_SE ? my $dN_VAR = $dN_SE*$dN_SE: 0;
                        $dS_SE ? my $dS_VAR = $dS_SE*$dS_SE: 0;
                        
                        # If there is no result (empty values of dN) then, skip this line
                        if($N || $S || $dS){
                        
                            push(@dN     ,$dN);
                            push(@dS     ,$dS);
                            push(@dN_VAR ,$dN_VAR);
                            push(@dS_VAR ,$dS_VAR);
                            push(@N      ,$N);
                            push(@S      ,$S);
                            push(@omega  ,$omega);
                            push(@kappa  ,$kappa);
                            push(@t      ,$t);
                            $out_str.= "$N\t$S\t$dN\t$dN_SE\t$dS\t$dS_SE\t$omega\t$kappa\t$t\n";
                        
                        }else{
                            
                            # We must discuss this awfull solution!!
                            
                            push(@dN     ,0);
                            push(@dS     ,0);
                            push(@dN_VAR ,0);
                            push(@dS_VAR ,0);
                            push(@N      ,0);
                            push(@S      ,0);
                            push(@omega  ,0);
                            push(@kappa  ,0);
                            push(@t      ,0);
                            $out_str.= "??\t??\t??\t??\t??\t??\t??\t??\t??\n";
                            
                        }
                    }
            
                }
            
            }
     
        }
        
    $out_str   .= "----------------------------------------------------------------------\n";
        
    my %hash = (
                yang_table => $out_str,
                N          => \@N,
                S          => \@S,
                dN         => \@dN,
                dS         => \@dS,
                dN_VAR     => \@dN_VAR,
                dS_VAR     => \@dS_VAR,
                omega      => \@omega,
                kappa      => \@kappa,
                t          => \@t
                );
    
    return (%hash);

}1;
