#!/usr/bin/env perl

use strict ;
use warnings ;
use lib "/usr/molgentools/lib";
use anne_files ;


# ---------------------------------------------------------- parameters -------------------------------------------------------------------------

my $tmpdrive = '/tmpdrive/bagel4wrapper/';
my $usage = "sudo /data/bagel4/bagel4_clean_tmpfolders.pl";


# get the folder sizes

my $command = "du --max-depth=1 /tmpdrive/bagel4wrapper/ >folder_sizes.txt" ;
system($command) ;
my @lines = anne_files::read_lines("folder_sizes.txt") ;
my $count = 0;
my @result ;
foreach my $line (@lines) {
	if ($line =~ m/(\d+)\s+(.*)/) {
		if ($1 <100) {
			$count++;
			push @result, "sudo -u www-data rm -r $2\n";
		}
	}
}
anne_files::write_lines("remove_bagel4_folders.sh", @result) ;
print "$count folders <100k\n" ;
my $totalFolders = scalar @lines ;
print "$totalFolders total folders\n";
print "sudo -u www-data ./remove_bagel4_folders.sh\n";

# sudo -u www-data rm -r /tmpdrive/bagel4wrapper/193.219.81.2gcsp2g8e3gnfquo96ai53ucc9783
# sudo -u www-data ls /tmpdrive/bagel4wrapper/193.219.81.2gcsp2g8e3gnfquo96ai53ucc9783

