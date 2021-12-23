#!/usr/bin/env perl

# routine to update the Bacteriocin databases from TABLE to FASTA

use strict ;
use warnings ;
use lib "/usr/molgentools/lib";
use anne_files ;
use lib "/data/bagel4/lib" ;
use bagel4_functions ;

# ---------------------------------------------------------- main -------------------------------------------------------------------------
my %conf = bagel4_functions::read_conf("../bagel4.conf") ;

# convert the table to a combined hmm file
make_hmm_db("../$conf{bacteriocin_hmm_table}", "bacteriocin_hmm") ;	


# ---------------------------------------------------------- functions -------------------------------------------------------------------------

sub make_hmm_db {
	my $infile = shift ;
	my $outfile = shift ;
	my %table = anne_files::read_table_to_hash($infile) ;
	my @HMMs ;
	foreach my $name (sort keys %table) {  
		print "Updating $name\n";
		push @HMMs, anne_files::read_lines("$name.hmm") ;
	}
	anne_files::write_lines($outfile,@HMMs) ;

}



