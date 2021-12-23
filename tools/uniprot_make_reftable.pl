#!/usr/bin/env perl

#  BAGEL3 program 
# 	Anne de Jong
#	University of Groningen
#	the Netherlands
#	anne.de.jong@rug.nl
#
#  2012-november
#
#  Makes a reference table for uniprot IDs
# 

use strict ;
use warnings ;

my $uniprotfile = "/var/uniref50/uniref50.fasta";
open FILE, "<$uniprotfile" or die "could not find $uniprotfile\n";
my @lines = <FILE>;	
close FILE ;
print "File loaded, start parsing...\n";
my $outputfile = "uniReference.db";
# example of header: >UniRef50_P39090 Alpha-1B-glycoprotein (Fragment) n=2 Tax=cellular organisms RepID=A1BG_EQUAS
open FILE, ">$outputfile" or die "could write to $outputfile\n";
# print header
print FILE "uniref_id\tdescription\n";
foreach my $line (@lines) {
	if ($line =~ m/\>(.*)\s(.*)\sn\=/g) {
		print FILE $1."\t".$2."\n";
	}	
}
print "Done\n";

