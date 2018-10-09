#!/usr/bin/perl

# construct input files for happy from raw data in CHR1

# Get all the genotypes.

$useful  =  "marker.order";
$output  =   "real_locus.data";
$alleles =  "real_locus.alleles";
$strains =  "strains";

## Deal with options
use Getopt::Std ;

getopts("ha:m:p:o:s:");        
if ( $opt_h ) { &usage } 

if ($opt_m) { 
$useful = $opt_m;
}

else {
	$opt_m = "marker.order";
	}

if ($opt_p) {
    $phenofile=$opt_p;
}
else { 
    $phenofile = "emo";
}
if ($opt_o) {
    $output=$opt_o;
}
if ($opt_s) {
    $strains=$opt_s;
}
if ($opt_a) {
	$alleles = $opt_a;
	}  


# Just use the useful markers

open(USEFUL, $useful) || die "Could not open $useful\n";
$n = 0;
while(<USEFUL>) {
    chomp;
    if ( /^\s*(\S+)\s+(\S+)/) {
	$useful{$1} = 1;
	$position{$1} = $2;
	$marker_order[$n++] = $1;
	#printf " Marker $marker_order[$n-1] position $position{$1}\n";
    }
}

close(USEFUL);
# Raw genotypes

foreach $marker ( keys %useful ) {
    my $mfile = "$marker";
    my $n = 0;
    if ( open( MARKER, $mfile ) ) {
	warn "$mfile\n";
	while(<MARKER>) {
	    chomp;
	    my ($mouse, $g1, $g2) = split;
	    $g1{$mouse,$marker} = $g1;
	    $g2{$mouse,$marker} = $g2;
	    $count{$marker,$g1}++;
	    $count{$marker,$g2}++;
	    $mouse{$mouse}++;
	    $order[$n++] = $mouse;
	}
	close(MARKER);
    }
    else {
	warn "Could not open $mfile\n";
	delete $useful{$marker};
	delete $position{$marker};
    }
}

# The phenotypes 

$phenotypes = "$phenofile";
open(PHENOTYPES,$phenotypes ) || die "Could not open $phenotypes\n";
$n = 0;
while(<PHENOTYPES>) {
    chomp;
    $phenotype{$order[$n++]} = $_;
}

# The strains information



open(STRAINS, $strains ) || die "Could not open $strains\n";

$_ = <STRAINS>;
chomp;
my @strains = split( /\s+/ );
shift @strains;
#shift @strains;

$m=0;
while(<STRAINS>) {
    chomp;
    my @data = split( /\s+/ );
    my $marker = $data[0];
#	print "@data\n";
    if ( $useful{$marker} ) {
	shift @data;
#	shift @data;

	if ( $#data == $#strains ) {
#	    $order[$m++] = $marker;
	    #print "strain $marker\n";
	    $n = 0;
	    foreach $allele ( @data ) {
		push @{$allele_list{$marker}},  $allele if ( ! $redundancy{$marker,$allele} );
		$redundancy{$marker,$allele}++;
		$allele{$marker,$strains[$n++]} = $allele;
	    }
	    $useful{$marker}=2;
	}
	else {
	    warn "$marker: incomplete strain data $#data $#strains\n";
	    delete $useful{$marker};
	}
    }
}

foreach $marker ( keys %useful ) {
    delete $useful{$marker} if ( $useful{$marker} != 2 );
} 

# check for errors

foreach $marker ( @marker_order ) {
    if ( $useful{$marker} ) {
	my %error = ();
	foreach $mouse ( keys %mouse ) {
	    if ( ($allele = $g1{$mouse,$marker} ) && ( !exists($redundancy{$marker,$allele}))) {
		$error{$allele} ++;
	    }
	    if ( ($allele = $g2{$mouse,$marker} ) && ( !exists($redundancy{$marker,$allele})) ) {
		$error{$allele} ++;
	    }
	}
	foreach $allele ( keys %error ) {
	    warn "marker $marker error allele $allele $error{$allele}  $redundancy{$marker,$allele}\n" if ( $allele ne "ND" );
	}
	foreach $key ( keys %redundancy ) {
	    ($m,$a) = split ( $; , $key );
	    if ( $m eq $marker ) {
		if ( ! exists($count{$m,$a}) ) {
		    warn "marker $marker missing alleles of type $a\n";
		}
	    }
	}
    }
}


# write out the genotypes
open( GENOS, ">$output") || die "Could not open $output\n";

$markers = 0;
printf GENOS "# %-15s ", "";

foreach $marker ( @marker_order ) {
    if ( $useful{$marker} ) {
	$markers++;
	printf GENOS " %9s", $marker;
    }
}
printf GENOS "\n";

$mice = 0;
foreach $mouse ( keys %mouse ) {
    $mice++;
    printf GENOS "%-15s $phenotype{$mouse} ", $mouse;
    foreach $marker ( @marker_order ) {
	if ( $useful{$marker} ) {
	    if ( $g1{$mouse,$marker} ) {
		printf GENOS " %4s", $g1{$mouse,$marker};
	    }
	    else  {
		printf GENOS " %4s", "ND";
	    }
	    if ( $g2{$mouse,$marker} ) {
		printf GENOS " %4s", $g2{$mouse,$marker};
	    }
	    else  {
		printf GENOS " %4s ", "ND";
	    }
	}
    }
    printf GENOS "\n";
}

close(GENOS);

warn "mice: $mice markers: $markers\n";

# write out the strains

open(ALLELES, ">$alleles") || die "Could not open $out\n";

print ALLELES "markers " , $markers , " strains ", $#strains+1, "\n";
print ALLELES "strain_names ", join(" ", @strains ), "\n";
foreach $marker ( @marker_order ) {
    if ( $useful{$marker} ) {
	$len = $#{$allele_list{$marker}}+2;
	print ALLELES "marker $marker $len $position{$marker}\n";
	printf ALLELES "allele %4s ", "ND";
	foreach $strain ( @strains ) {
	    printf ALLELES " %.3f", 1/($#strains+1);
	}
	print ALLELES "\n";

	foreach $allele ( @{$allele_list{$marker}} ) {
	    printf ALLELES "allele %4s ", $allele;
	    foreach $strain ( @strains ) {
		if ( $allele{$marker,$strain} eq $allele ) {
		    printf ALLELES " %.3f", 1/$redundancy{$marker,$allele};
		}
		else {
		    printf ALLELES " %.3f", 0.0;
		}
	    }
	print ALLELES "\n";
        }
    }
}

close(ALLELES);

#------------------------------------------------------------------------------

sub usage {

die "
        qtlData.pl <options> 

       Looks for the following files:
             marker.order - two columns, marker and distance. Markers in order
		The names of the markers in this file are used to search in the 
		current directory for the corresponding genotype files        

	     strains      - allele sizes in bp for each marker. # lines are ignored.
			eg D1MIT29	234	230	230	228 etc
	
        OPTIONS:


                -h                      : usage
                
		-o			: output locus data file name (default: real_locus.data)
		-a			: output alleles file name    (default: real_locus.alleles)
		-s                      : input strain file name      (default: strains)
		-m			: input marker order file name(default: marker.order)   
		-p			: input phenotype file
" ;

}          
