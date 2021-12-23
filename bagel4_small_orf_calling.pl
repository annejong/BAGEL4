#!/usr/bin/env perl

#  BAGEL4 program 
# 	Anne de Jong and Auke van Heel
#	University of Groningen
#	the Netherlands
#	anne.de.jong@rug.nl
#
#  	samll ORF calling routine and filter for Lantibiotics properties 
#
# v 1.0 
# 


use strict ;
use warnings ;
use lib "/data/bagel4/lib";

use bagel4_functions ;
use lib "/data/molgentools/lib";
use anne_files ;
use anne_genomics ;

my %conf = bagel4_functions::read_conf('/data/bagel4/bagel4.conf') ;


# -----------------------------------------------------------------------------------  global vars ---------------------------------------------------------------------------
my $sessiondir = "/data/bagel4/test" ;
my $queryname ;
my $min = 60 ;
my $max = 180 ;
my $ORF_overlap = 30 ;  # allows 10bp overlap
my $usage = "bagel4_small_orf_calling.pl
				-s Sessiondir [default=$sessiondir]
				-queryname queryfile
				-min minimal aa [default=$min]
				-max maximal aa [default=$max]
				
		DEPENDS on: queryname.fna & queryname.predict
		
		small ORF calling routine 		
				
e.g.  /data/bagel4/bagel4_small_orf_calling.pl -s /tmp/BAGEL4WRAPPER/129.125.143.35q44h5a3a76ka7353q1dd0nhl45828 -queryname NC_0085331.AOI_1 -min 60 -max 240 
	
" ;

&parseparam() ;


# ---------------------------------------------------------- main -------------------------------------------------------------------------


# 1. Small orf calling
	my $DNA = uc(anne_genomics::get_dna_from_fasta("$sessiondir/$queryname.fna")) ;
	my $DNA_masked = mask_orf_regions($DNA) ;
	my $orf_count = 0 ;
	my %GeneTable ;
	&orf_calling($DNA_masked, '+');
	&orf_calling(anne_genomics::inverse_complement($DNA_masked), '-');
	
	my @results_fna ;
	my @results_table = "Queryname\tstart\tend\tstrand\trbs\tgene" ;
	#print "SMALL ORFs:\n" ;
	#print "\tQueryname\tstart\tend\tstrand\trbs\n";
	foreach my $key (sort {$GeneTable{$a}{start} <=> $GeneTable{$b}{start}} keys %GeneTable) {
		push @results_fna, ">$key strand=$GeneTable{$key}{strand} start=$GeneTable{$key}{start} end=$GeneTable{$key}{end} RBS=$GeneTable{$key}{rbs}\n$GeneTable{$key}{gene}\n"  ; 
		push @results_table, "$key\t$GeneTable{$key}{start}\t$GeneTable{$key}{end}\t$GeneTable{$key}{strand}\t$GeneTable{$key}{rbs}\t$GeneTable{$key}{gene}" ;
		my $protein = anne_genomics::translate($GeneTable{$key}{gene});
		my $dna = substr $GeneTable{$key}{gene}, 0, 30 ;
		my $dna_true = substr $DNA_masked, $GeneTable{$key}{start}, 30 ;
		my $truepos = index($DNA_masked, $GeneTable{$key}{gene});
		$dna_true = substr anne_genomics::inverse_complement($DNA_masked), $GeneTable{$key}{start}, 30 if ($GeneTable{$key}{strand} eq '-') ;
		#print "\t$key\t$GeneTable{$key}{start}\t$truepos\t$GeneTable{$key}{end}\t$GeneTable{$key}{strand}\t$GeneTable{$key}{rbs}\t$dna\.\.\.  $dna_true\n";
	}	
	
	anne_files::write_lines("$sessiondir/$queryname.sOFS.fna", @results_fna) ;
	anne_files::write_lines("$sessiondir/$queryname.sOFS.table", @results_table) ;
	


# ---------------------------------------------------------- functions -------------------------------------------------------------------------

sub mask_orf_regions {
	# mask all gene regions and the tailing parts
	my $dna = shift ;
	my @lines = anne_files::read_lines("$sessiondir/$queryname.predict");  # read the glimmer orf positions
	my $first_start = 0 ;
	my $last_end = 0;
	my $first_line = 1 ;
	foreach my $line(@lines) {
		my @items = split /\s+/, $line ;
		if (scalar @items>2) {
			if ($items[1] > $items[2]) { my $tmp = $items[1] ; $items[1] = $items[2]; $items[2] = $tmp; }  # swap if start/end for strand -
			my $start = $items[1] + $ORF_overlap ; # allow ORF overlap
			my $end   = $items[2] - $ORF_overlap ;
			if ($first_line) { $first_line = 0; $first_start = $start; }
			$last_end = $end;
			my $len = $end-$start ;
			my $N_ORF = 'N' x $len ;
			my $ORF = substr $dna, $start, $len ;
			$dna =~ s/$ORF/$N_ORF/ ;  # put N's at the ORF regions
		}	
	}
	
	# Do not allow ORF calling in the tails; remove the tails
	my $N_ORF = 'N' x $first_start ;
	my $ORF = substr $dna, 0, $first_start ;
	$dna =~ s/$ORF/$N_ORF/ ;  # put N's at the left tail
	my $len = length($dna) - $last_end  ;
	#print "$N_ORF\n$ORF\n$len-$last_end\n";
	if ($len>1) {
		$N_ORF = 'N' x $len ;
		$ORF = substr $dna, -$len ;
		$dna =~ s/$ORF/$N_ORF/ ;  # put N's at the right tail
		#print "$N_ORF\n$ORF\n$len-$last_end\n";
	}	
	return $dna ;
}


sub orf_calling {
	my $DNA = shift ;
	my $strand = shift ;
	my $startcodons = 'ATG|GTG|TTG';
	my $index = 1 ;
	my $DNAlen = length($DNA) ;
	my $start ;
	my $end ;
	my $lastgene ="";
	#	while ( (substr $DNA, $index) =~m/(..GG.{5,14})((ATG|GTG)(...){0,$max}?(TGA|TAG|TAA))/) {
	while ( (substr $DNA, $index) =~m/(.....{5,14})(($startcodons)((ATG|GTG|TTG|GCC|AGT|TGT|CGA|ATC|AAC|AGC|TAC|TCG|ACA|CTG|CCG|GCA|AAG|GTT|CAC|AGA|ACC|CCA|TGG|CGC|CTC|CAG|ACG|AAA|GTA|CTT|GGA|GTC|TGC|TCA|ATT|TAT|AAT|ACT|CAA|GAC|GGT|TCC|TTT|AGG|CGT|CGG|CAT|ATA|CCC|GGG|GAG|TTA|CTA|GAT|TCT|TTC|GCG|GGC|GAA|GCT|CCT)){$min,$max}(TGA|TAG|TAA))/) {
		my $RBS = $1 ;
		my $gene = $2 ;
		my $currentIndex = $index ;
		my $zeroIndex = $-[0] ;
		my $genestart = $index + length($RBS) + $-[0] ;
		$index += $-[0] + 3*$conf{smallorf_min_prot_length} - 3 ;  # goto the next position and check again for a gene
		#my $genestart = $index - 3*$conf{smallorf_min_prot_length} - 3 ;  # the gene start position
		my $DNAstart = substr $DNA, $genestart, 10 ;
		$RBS =~ s/N//g ;  # remove the N's from the RBS after the $_ values because s/ will generate a new $_
		if ($gene ne $lastgene) {
			$orf_count++ ;
			my $key = "sORF_$orf_count" ;
			#print "$key\t$strand\t$currentIndex\t$index\t$zeroIndex\t$genestart\t$DNAstart\t".anne_genomics::translate($gene)."\n";
			$GeneTable{$key}{gene} = $gene ;
			$GeneTable{$key}{rbs} = $RBS ;
			if ($strand eq '+') {
				$GeneTable{$key}{start}	= $genestart + 1 ;
				$GeneTable{$key}{end}	= $GeneTable{$key}{start} + length($gene) -1;
				$GeneTable{$key}{strand} = $strand ;
			} elsif ($strand eq '-') {
				$GeneTable{$key}{end} 	= $DNAlen - $genestart ;
				$GeneTable{$key}{start}	= $GeneTable{$key}{end} - length($gene) + 1 ;
				$GeneTable{$key}{strand} = $strand ;
			}	
		}
		$lastgene = $gene ;
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
        die $usage if ($var eq '-h' or $var eq '--help') ;
		$sessiondir		= shift(@arg) if($var eq '-s') ;
		$queryname		= shift(@arg) if($var eq '-queryname') ;
		$min		= shift(@arg) if($var eq '-min') ;
		$max		= shift(@arg) if($var eq '-max') ;
    }
	die $usage if (!$queryname) ;
}

