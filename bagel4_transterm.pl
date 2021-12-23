#!/usr/bin/env perl


# Transcription terminator prediction routine adapted from genome2D

# transtermHP


use strict ;
use warnings ;
use lib "/usr/bagel4/lib" ;
use bagel4_functions ;
use lib "/usr/molgentools/lib";
use anne_files ;
use anne_genomics ;
use File::Basename;


# ---------------------------------------------------------- parameters -------------------------------------------------------------------------

my $sessiondir = '/usr/bagel4/test';
my $program_dir = dirname($0) ;
my $fna_file = 'AOI.fna';
my $usage = "option:
		-s		sessiondir [default=$sessiondir]
		-prefix	prefix of the input : .fna .ptt .coords 
				
		Transcription terminator prediction routine adapted from genome2D
		output is tab delimtied table: position strand sequence
		NOTE: Fasta file must have the same header name as the prefix and must be in the \$sessiondir
		
				
		e.g.  /usr/bagel4/bagel4_transterm.pl -s /usr/bagel4/test -dna NC_0085331.AOI.fna
" ;


&parseparam() ;
my %conf = bagel4_functions::read_conf("$program_dir/bagel4.conf") ;

# ---------------------------------------------------------- main -------------------------------------------------------------------------

my $pvalue = 0.00005 ;

my %fasta = anne_genomics::fasta2hash($fna_file) ;  # read the multi fasta file
foreach my $key (keys %fasta) {
	my $command = "$conf{transtermHP}/transterm --bag-output $sessiondir/$key.transterm.results -p $conf{transtermHP}/expterm.dat $sessiondir/$key.fna $sessiondir/$key.ptt > $sessiondir/transterm.log";
	print "$command\n" ;
	system ($command) ;
	anne_files::write_lines("$sessiondir/$key.transterm", parse_transterm_results("$sessiondir/$key.transterm.results") );	
}

# ---------------------------------------------------------- functions -------------------------------------------------------------------------



sub parse_transterm_results {
# e.g. result: NC_003098.1 Streptococcus pneumoniae R6 chromosome, complete genome,SigmaN16.pfm,110292,-,10.6753806385,ATATTAAAGAGACTCAAAAAGTTGTCAA,
	# write the position and strand to file
	my $filename = shift ;
	my @results ;
	my @lines = anne_files::read_lines($filename) ;
	foreach my $line (@lines) {
		if ($line !~ m/NONE/g) {
			$line =~ s/^\s+// ; # remove the starting spaces
			my @items = split /\s+/, $line ;
			if (scalar @items > 4) { 
				push @results, "$items[1]\t$items[4]\t$line" ;
				#print "$items[1]\t$items[4]\t$line\n" ;
			}
		}	
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


