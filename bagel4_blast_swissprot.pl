#!/usr/bin/env perl

# Anne de Jong
# Uniprot blast routine for BAGEL4, Nov 2017

use strict ;
use warnings ;
use lib "/data/molgentools/lib";
use anne_files ;
use File::Basename;

my $sessiondir = '/data/bagel4/swissprot';
my $program_dir = dirname($0) ;
my $query ;
my $db = '/data/bagel4/swissprot/archea_bacteria' ;
my $cpu = 4 ;
my $evalue = 0.00001 ;
my $blastresultfile = 'blast.results';
my $outfile = 'blast_uniprot.table' ;



my $usage = "options:
		-s	sessiondir [default=$sessiondir]
		-query	multiple entries fasta
		-db	protein database (formatted)
		-cpu	number of cpu's used [default=$cpu]
		-evalue 
		-o	output table name [default=$outfile]
				
		this routine will perform blast all proteins in the query to the database $db, see README in the /swissprot folder for details
		blastp is used to indentify the proteins
		descriptions of the full db fasta headers are added to the results
				
		e.g.  time /data/bagel4/blast_swissprot.pl -s /data/bagel4/swissprot -db archea_bacteria -query /data/bagel4/swissprot/example_opp.fasta -name my_query -cpu 4 
" ;

parseparam();
# --------------------------------------------- main -------------------------------------------------------------------------------------

blastp();

my @results = add_annotation() ;

anne_files::write_lines("$outfile", @results) ;
print "Results written to $outfile\n";

# --------------------------------------------- functions -------------------------------------------------------------------------------

sub add_annotation {
	my %descriptions = anne_files::read_table_to_hash("$program_dir/swissprot/bagel4_swissport_descriptions_pfams.table") ;
	my %blast_results = anne_files::read_table_to_hash("$sessiondir/blast_results_file") ;
	my @result = "Query\tUniProt\tPercent\tEvalue\tFunction\tTax\tPFAMS\turl";
	foreach my $key (sort keys %blast_results) {
		my $function = '';	$function = $descriptions{$blast_results{$key}{Subject}}{Function} if (defined($descriptions{$blast_results{$key}{Subject}}{Function}));
		my $Tax = '';		$Tax =      $descriptions{$blast_results{$key}{Subject}}{Tax} if (defined($descriptions{$blast_results{$key}{Subject}}{Tax})) ;
		my $pfams = ''; 	$pfams =    $descriptions{$blast_results{$key}{Subject}}{PFAMS} if (defined($descriptions{$blast_results{$key}{Subject}}{PFAMS})) ;
		#$blast_results{$key}{Subject} =~ s/UniRef90_// ;
		my $url = "<a href=http://www.uniprot.org/uniprot/$blast_results{$key}{Subject}>$blast_results{$key}{Subject}</a>";
		push @result, "$key\t$blast_results{$key}{Subject}\t$blast_results{$key}{percent}\t$blast_results{$key}{evalue}\t$function\t$Tax\t$pfams\t$url";
		#print "$key\t$blast_results{$key}{Subject}\t$blast_results{$key}{percent}\t$blast_results{$key}{evalue}\t$function\t$Tax\t$pfams\t$url\n";
	}
	return @result ;
}


sub blastp {
	my $command = "blastp -num_threads $cpu -query $query -db $db -outfmt 6 -max_target_seqs 1 -evalue $evalue -out $sessiondir/blast_results_file" ;
	print "===>$command\n"; 
	system($command) ;
	&anne_files::add_header("$sessiondir/blast_results_file", "Query\tSubject\tpercent\talign_len\tmistmatches\tgap\tq_start\tq_end\ts_start\ts_end\tevalue\tbit_score");
}



sub parseparam {
    my $var ;
    my @arg = @ARGV ;
    while(@arg) {
        $var = shift(@arg) ;
        die $usage if ($var eq '-h' or $var eq '--help') ;
		$sessiondir	= shift(@arg) if($var eq '-s') ;
		$query		= shift(@arg) if($var eq '-query') ;
		$db			= shift(@arg) if($var eq '-db') ;
		$evalue		= shift(@arg) if($var eq '-evalue') ;
		$cpu		= shift(@arg) if($var eq '-cpu') ;
		$outfile	= shift(@arg) if($var eq '-out') ;
    }
    die $usage if (!$query) ;
}