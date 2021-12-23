#!/usr/bin/env perl

# routine to update the Bacteriocin databases from TABLE to FASTA

use strict ;
use warnings ;
use lib "/usr/molgentools/lib";
use anne_files ;


# ---------------------------------------------------------- main -------------------------------------------------------------------------

# convert the tables
#table_2_fasta('bagel4_class1bacteriocin_db.txt','bagel4_class1bacteriocin_db.fasta');
table_2_fasta('bagel4_class2bacteriocin_db.txt','bagel4_class2bacteriocin_db.fasta');
#table_2_fasta('bagel4_class3bacteriocin_db.txt','bagel4_class3bacteriocin_db.fasta');

# format fasta for BLAST
my $command = 'formatdb -p T -o T -i bagel4_class1bacteriocin_db.fasta' ; system($command) ;
   $command = 'formatdb -p T -o T -i bagel4_class2bacteriocin_db.fasta' ; system($command) ;
   $command = 'formatdb -p T -o T -i bagel4_class3bacteriocin_db.fasta' ; system($command) ;

# ---------------------------------------------------------- functions -------------------------------------------------------------------------

sub table_2_fasta {
	my ($infile, $outfile) = @_ ;
	my %table = anne_files::read_table_to_hash($infile) ;
	my @lines ;
	foreach my $ID (sort keys %table) {
		if (defined($table{$ID}{Sequence}) and length($table{$ID}{Sequence})>2)	 {
			push @lines, ">$ID;$table{$ID}{NAME}";
			push @lines, $table{$ID}{Sequence};
		} else {
			print "No Sequence found for $ID\t$table{$ID}{Sequence}\n";
		}	
	}
	anne_files::write_lines($outfile, @lines);
}	



