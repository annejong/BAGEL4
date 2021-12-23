#!/usr/bin/env perl

# Anne de Jong
# 
# Written for BAGEL4
# Raw orf calling routine

use strict ;
use warnings ;
use lib "/usr/bagel3/lib";
use genomics ;

# -----------------------------------------------------  parameters  ---------------------------------------
my $sessiondir = './';
my $inputfile ;
my $min = 120;
my $max = 1500 ;
my $outputfile = "results.txt";
my $RBSregex = '..GG.{5,14}';
my $startcodons = 'ATG|GTG|TTG';
my $usage = "./orf_calling
				-s Sessiondir [default=current folder]
				-i fasta DNA file
				-min minimal protein size [default=$min]
				-max maximal protein size [default=$max]
				-rbs [default=..GG.{5,14}]
				-startcodons [default=ATG|GTG|TTG]
				-o output file [default=results.txt] two result files: DNA_<filename>  and PROT_<filename> \n
		e.g.  /usr/bagel3/orf_calling.pl -i NC_009004.fna 		
" ;

&parseparam() ;
print "ORF calling on $inputfile with min=$min aa, max=$max aa\n";
my $genome = &read_DNA() ;
	open (TMP, '>',  "$sessiondir/QueryDNA.fna") ;
	print TMP "$genome\n" ;
	close TMP ;

print "genome size= ".length($genome)." bp\n";
my $orf_count = 0 ;
my $prot_count = 0 ;
&find_orfs_both_strands();



# ---------------------------------------------------------------------------------------- main ----------------------------------------------------------------------------------------

sub find_orfs_both_strands {
	open (DNA, '>',  "$sessiondir/DNA_$outputfile") or die ("unable to write $sessiondir/DNA_$outputfile\n") ;
	open (PROT, '>',  "$sessiondir/PROT_$outputfile") or die ("unable to write $sessiondir/PROT_$outputfile\n") ;
	print "Scanning upper strand\n";
	&find_orfs_with_overlap_fast($genome, 'strand+');
	print "Scanning lower strand\n";
	&find_orfs_with_overlap_fast(&genomics::inverse_complement($genome), 'strand-');
	print "$orf_count genes\n";
	print "$prot_count proteins fits criteria\n";
	close DNA ;
	close PROT ;
	if ($orf_count == 0) {
		unlink ("$sessiondir/DNA_$outputfile") ;
		unlink ("$sessiondir/PROT_$outputfile") ;
	}	
}

sub find_orfs_with_overlap_fast {
	my $DNA = shift ;
	my $ori = shift ;
	my $index = 1 ;
	my $DNAlen = length($DNA) ;
	my $start ;
	my $end ;
	my $lastgene ="";
	my $mingene= $min*2;
	# remove the /g option to get the first hit only
#	while ( (substr $DNA, $index) =~m/(..GG.{5,14})((ATG|GTG)(...){0,$max}?(TGA|TAG|TAA))/) {
	while ( (substr $DNA, $index) =~m/(..GG.{5,14})(($startcodons)((ATG|GTG|TTG|GCC|AGT|TGT|CGA|ATC|AAC|AGC|TAC|TCG|ACA|CTG|CCG|GCA|AAG|GTT|CAC|AGA|ACC|CCA|TGG|CGC|CTC|CAG|ACG|AAA|GTA|CTT|GGA|GTC|TGC|TCA|ATT|TAT|AAT|ACT|CAA|GAC|GGT|TCC|TTT|AGG|CGT|CGG|CAT|ATA|CCC|GGG|GAG|TTA|CTA|GAT|TCT|TTC|GCG|GGC|GAA|GCT|CCT)){$min,$max}?(TGA|TAG|TAA))/) {
		my $RBS = $1 ;
		my $gene = $2 ;
		$index += $-[0] + 1 ;  # goto the next base and check again for a gene
		#if (length($gene)>=$mingene) 
		{
			my $genestart = $index + length($RBS) -1;  # the gene start position
			if ($gene ne $lastgene) {
				$orf_count++ ;
				if ($ori eq 'strand+') {
					$start = $genestart + 1 ;
					$end   = $start + length($gene) -1;
				} elsif ($ori eq 'strand-') {
					$end = ($DNAlen-$genestart ) ;
					$start  = $end - length($gene) + 1 ;
				}	
				print DNA ">ORF_$orf_count $ori start=$start RBS=$1\n$gene\n"  ;
				my $prot = &genomics::translate($gene) ;
				$prot_count++;
				$prot =~ s/\-//g ;
				print PROT ">ORF_$orf_count $ori LEN=".length($prot)."aa start=$start end=$end RBS=$RBS GENE=gene min=$min max=$max ".length($gene)."\n" ;
				print PROT "$prot\n";
			}
			$lastgene = $gene ;
		}	
	}
}



sub chk_protein {
	my $prot = shift ;
	my $chk = 0 ;
	if ( (substr $prot, 1, 15) !~ m/.*C.*/g ) {  # should not contain a C in the first 15 aa
		my $mature = substr $prot, 15  ;
		if ( $mature=~ m/[C]/g and $mature=~ m/[TS]/g ) {  # should contain a C or T/S in the leader
			$chk = 1;
		}	
	}
	return $chk ;
}

sub parseparam {
    my $var ;
    my @arg = @ARGV ;
    while(@arg) {
        $var = shift(@arg) ;
        $sessiondir = shift(@arg) if($var eq '-s') ;
        $inputfile 	= shift(@arg) if($var eq '-i') ;
		$min		= shift(@arg) if($var eq '-min') ;
		$max		= shift(@arg) if($var eq '-max') ; 
		$RBSregex	= shift(@arg) if($var eq '-rbs') ; 
		$startcodons= shift(@arg) if($var eq '-startcodons') ; 
        $outputfile = shift(@arg) if($var eq '-o') ;
    }
    die "No filename found\n$usage" if (!$inputfile) ;

}

sub read_DNA {
	# load the inputfile
	open (FILE,'<', "$inputfile") or die ("$inputfile does not exists") ;
	my(@lines) = <FILE>;					# Get the tablefile into a array
	chomp @lines;
	close FILE ;
	my $header ;
	my $DNA ;
	foreach my $line (@lines)  { 
		$line =~ s/(\n|\r|\x0d)//g;  # remove DOS/windows line breaks etc
		if ($line =~m/^\>/) {
			$header = $line ;
		} else {
			$DNA .= uc($line) ;
		}
	}
	return $DNA ;
}

