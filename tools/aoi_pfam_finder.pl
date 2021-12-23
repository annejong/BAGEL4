#!/usr/bin/env perl
use strict ;
use warnings ;
use lib "/usr/molgentools/genomics";
use genomics ;




my %condition = genomics::read_table_to_hash("hmm_db/bagel3.primary_pfam_rules.db");
genomics::print_hash_of_hashes(%condition);
my $AOIname = 'my_name' ;
$condition{$AOIname}{rule} = "(siger|anne) and (auke|heel)" ;
$condition{$AOIname}{AOIsize} = 800 ;
$AOIname = 'PF_test' ;
$condition{$AOIname}{rule} = "(PF04737|anne) and (PF05147)" ;
$condition{$AOIname}{AOIsize} = 10000 ;


# ---------- implementation -----------------

my @pfams = get_unique_pfams() ; # get a list of unique pfams from the rules of the condition table

exit ;

my %hmmresults = parse_hmmsearchresults("my_session/all.hmmprimary.result"); # key is the pfam name with a sequence number e.g. PF1234.1_2
genomics::print_hash_of_hashes(%hmmresults) ;
$hmmresults{anne}{position} = 1000 ;
$hmmresults{'auke.1'}{position} = 2000 ;
$hmmresults{'auke.2'}{position} = 2600 ;
$hmmresults{siger}{position} = 3000 ;



my %pfamAOIs = locate_AOI() ; # key is the pfam name on which the AOI is called
genomics::print_hash_of_hashes(%pfamAOIs) ;

sub get_unique_pfams {
	my %pfams ;
	foreach my $key (sort keys %condition) {
		$condition{$key}{rule} =~ s/ and / AND /g ;
		my @items = split ( / AND /, $condition{$key}{rule} ) ; 
		foreach my $item (@items) {
		print "==>$item\n" ;
			$item =~ s/(\(|\))//g ; # remove brackets
			my @sub_items = split ( /\|/, $item ) ; 
			foreach my $tmp (@sub_items) { $pfams{$tmp} = 1 ; }
		}
	}
	foreach my $key (sort keys %pfams) {
		print "$key\n";
	}
}

sub parse_hmmsearchresults {
	# parse the results from one or multiple hmmsearch3 result files
	my $filename = shift ;
	my %result ;
	my @lines = genomics::read_lines($filename) ;
	chomp @lines ;
	my $count = 0 ;
	my $last_start = '';
	my $last_end = '';
	
	foreach my $line (@lines) {
		my @items = split /\s+/, $line ;
		if ($line =~ m/.*strand(.).*start=(\d+).*end=(\d+)/) {	
			if ($last_start ne $2 and $last_end ne $3) {
				$count++;
				my $pfam = $items[3]."_$count" ; 
				$result{$pfam}{strand} = $1 ;
				$result{$pfam}{start} = $2 ;
				$result{$pfam}{end} = $3 ;
				$result{$pfam}{position} = ( $result{$pfam}{start} + $result{$pfam}{end} ) / 2 ;
				#print "PFAM=$pfam\tstrand$1\tstart=$2\tend=$3\n";
			}	
			$last_start = $2 ;
			$last_end = $3 ;
		}
	}
	return %result ;
}
 

sub locate_AOI {
    my %results ;
    foreach my $AOIname (sort keys %condition) {
		# push all the pfams of the hmm results into one teststring
		my %teststrings ;
		foreach my $key (keys %hmmresults) {
			my $tmp ;
			foreach my $key2 (keys %hmmresults) {
			$tmp .= $key2.' ' if (abs($hmmresults{$key}{position} - $hmmresults{$key2}{position}) <= $condition{$AOIname}{AOIsize}) ;
			}	
			$teststrings{$key}{position} = $hmmresults{$key}{position}  ;
			$teststrings{$key}{teststring} = $tmp ;
		}

		# check all the teststring with all the rules; the rules contain or as an | for the regular expression, and the AND will be tested iteratly by an if statement 
		my @rules = split ( (' and '|' AND '), $condition{$AOIname}{AOIrule} ) ;
		foreach my $key (sort keys %teststrings) {
			my $count = 0 ;
			#print "checking $key:\n" ;
			foreach my $rule (@rules) {
				$count++ if ($teststrings{$key}{teststring} =~ m/$rule/) ;
			}
			if (scalar @rules == $count) {
				print "AOI found on position $teststrings{$key}{position} of $key on the basis of '$teststrings{$key}{teststring}' with rules '$condition{$AOIname}{AOIrule}'\n"; 
				$results{$key}{AOIname} = $AOIname ;
				$results{$key}{AOIsize} = $condition{$AOIname}{AOIsize} ;
				$results{$key}{position} = $teststrings{$key}{position} ;
			} 
		}	
    }
    return %results ;
}


    