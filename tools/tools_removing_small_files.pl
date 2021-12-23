#!/usr/bin/env perl

# remove small seq files

use strict ;
use warnings ;
use lib "/usr/molgentools/genomics";
use genomics ;

# ls -s > /tmp/files

my $cutoffsize = 5000 ; # approximately the number of basis

my @lines = genomics::read_lines('/tmp/files') ;
my $count = 0 ;
my $basecount = 0;
foreach my $line (@lines) {
	if ($line =~ m/\s+(\d+)\s(.*)/) {
		my $file = "/var/Bacteria_DRAFT.fna/$2" ;
		my $filesize = -s $file;
		$basecount = $basecount + $filesize ;
		if ($filesize < $cutoffsize) {
			$count ++ ;
			print "Removing $file with size $filesize [$count]\n";
			unlink $file;
		}
	}
}
my $filecount = scalar @lines ;
print "Number of files removed:$count out of $filecount\n";
print "Number of basis left: $basecount\n";