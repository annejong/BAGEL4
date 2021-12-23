#!/usr/bin/env perl

# Tools to download the latest RefSeq summary and make a selectionlist for the webserver


use strict ;
use warnings ; 
use Data::Dumper ;
use lib "/usr/molgentools/lib";
use anne_files ;
use File::Basename;


# ---------------------------------------------------------- parameters -------------------------------------------------------------------------

my $program_dir = dirname($0) ;


# ------------------------------------------------- parameters ------------------------------------------------------------------------
my $refseq_summary_file = "$program_dir/00.latest_assembly_summary.txt" ;	# the latest summary file of NCBI
my $webfolder = "/data/www/bagel4";

# 1. Download latest summary
	print "Download latest summary, contains all genome information\n" ;
	my $tmp = "wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt -O $refseq_summary_file";
	print "$tmp\n";
	#system($tmp) ;
	
# 2. Read the table
	my @lines = anne_files::read_lines($refseq_summary_file);
	
	my @my_list ;
	foreach my $line (@lines) {
		if ($line !~ m/^#/) {
			my @items = split /\t/, $line ;
			$items[8] =~ s/strain=//;
			push @my_list, "$items[0]=$items[7] $items[8]" ;
		}
	}
	print 'Number of Genomes: '.(scalar @my_list)."\n";
	
	
# 3. write list to webserver folder	
	anne_files::write_lines("$webfolder/RefSeq_WGS_list.txt", @my_list) ;
	
	