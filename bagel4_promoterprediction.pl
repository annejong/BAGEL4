#!/usr/bin/env perl


# Promoter prediction routine adapted from genome2D

# depends on MOODS 


use strict ;
use warnings ;
use lib "/data/bagel4/lib" ;
use bagel4_functions ;
use lib "/data/molgentools/lib";
use anne_files ;
use anne_genomics ;
use File::Basename;


# ---------------------------------------------------------- parameters -------------------------------------------------------------------------

my $sessiondir = '/tmp/BAGEL4WRAPPER/test';
my $program_dir = dirname($0) ;
my $fna_file  ;
my $outfile = 'promoterscan.table' ;
my $usage = "option:
		-s		sessiondir [default=$sessiondir]
		-dna	multi fasta file e.g. AOI.fna 
		-o		results file [deafult = $outfile]
				
	simple Promoter prediction routine adapted from genome2D
	output is tab delimtied table: position strand sequence
				
	e.g.  /data/bagel4/bagel4_promoterprediction.pl -s /tmp/BAGEL4WRAPPER/test -dna NC_0085331.AOI.fna -queryname NC_0085331	
" ;
my %conf = bagel4_functions::read_conf("$program_dir/bagel4.conf") ;

&parseparam() ;

# ---------------------------------------------------------- main -------------------------------------------------------------------------

my %fasta = anne_genomics::fasta2hash($fna_file) ;  # read the multi fasta file

mask_orf_regions() ;  # mask the ORF regions (by N's) 

my $pvalue = 0.00005 ;
foreach my $key (keys %fasta) {
	anne_files::write_string("$sessiondir/$key.masked.fna", ">$key\n$fasta{$key}\n");
	my $command = "$conf{MOODS}/scripts/moods_dna.py -p $pvalue -m $program_dir/promoterprediction/SigmaN*.pfm  -s $sessiondir/$key.masked.fna -o $sessiondir/$key.MOODS" ;
	print "$command\n";
	system ($command) ;
	anne_files::write_lines("$sessiondir/$key.promoters", parse_MOODS_results("$sessiondir/$key.MOODS")); 
}	
		

# ---------------------------------------------------------- functions -------------------------------------------------------------------------


sub mask_orf_regions {
 	foreach my $key (keys %fasta) {
		my @lines = anne_files::read_lines("$sessiondir/$key.predict");
		foreach my $line(@lines) {
			my @items = split /\s+/, $line ;
			if (scalar @items>2) {
				if ($items[1] > $items[2]) { my $tmp = $items[1] ; $items[1] = $items[2]; $items[2] = $tmp; }
				my $len = $items[2]-$items[1]-20 ;
				my $N_ORF = 'N' x $len ;
				my $ORF = substr $fasta{$key}, $items[1]+10, $len-10 ;
				$fasta{$key} =~ s/$ORF/$N_ORF/ ;  # put N's at the ORF regions
			}	
		}
	}	
}

sub parse_MOODS_results {
# e.g. result: NC_003098.1 Streptococcus pneumoniae R6 chromosome, complete genome,SigmaN16.pfm,110292,-,10.6753806385,ATATTAAAGAGACTCAAAAAGTTGTCAA,
	# write the position and strand to file
	my $filename = shift ;
	my @results ;
	my @lines = anne_files::read_lines($filename) ;
	foreach my $line (@lines) {
		my @items = split ",", $line ;
		if (scalar @items > 4) { 
			$items[5] = anne_genomics::inverse_complement($items[5]) if ($items[3] eq '-') ;
			push @results, "$items[2]\t$items[3]\t$items[5]" ; }
	}
	return @results ;
}



sub parseparam {
    my $var ;
    my @arg = @ARGV ;
    while(@arg) {
        $var = shift(@arg) ;
        die $usage if ($var eq '-h' or $var eq '--help') ;
		$sessiondir		= shift(@arg) if($var eq '-s') ;
		$fna_file		= shift(@arg) if($var eq '-dna') ;
		
    }
	die $usage if (!$fna_file) ;
}


