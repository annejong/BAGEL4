#!/usr/bin/env perl

use strict ;
use warnings ;
use lib "/data/bagel4/lib" ;
use bagel4_functions ;
use lib "/data/molgentools/lib";
use anne_files ;
use File::Basename;


# ---------------------------------------------------------- parameters -------------------------------------------------------------------------

my $sessiondir = '/data/bagel4/test';
my $program_dir = dirname($0) ;
my $AOI_identification_rules_table = "$program_dir/tables/AOI_identification_rules_v2.txt" ;
#my $AOI_identification_rules_table = "$program_dir/tables/AOI_identification_rules_v2_testing.txt" ;
my $outputfile = 'hmmsearch.results';
my $queryname = 'queryname';
my $fna_file ;
my $debug = 1 ;

# default settings: 
	my %conf = bagel4_functions::read_conf("$program_dir/bagel4.conf") ;
	my $contextsize = $conf{contextsize};
	my $domT = $conf{domT} ; # PFAM score cutoff --domT
	$domT = 10 ;
	my $cpu = $conf{cpu} ;
	my $distance = $conf{contextsize} ;

	
my $usage = "options:
	-s Sessiondir [default=$sessiondir]
	-dna dna input file in FASTA format
	-cpu number of cpu's used [default=$cpu]
	-queryname name of the query
	-size contextsize [default=$contextsize]
	-o Results file [default=$outputfile]
	
	this routine uses the identification rules from table AOI_identification_rules.txt [default=$AOI_identification_rules_table]
				
		e.g.  /data/bagel4/bagel4_AOI_identification_hmm_rules.pl -dna /data/genomes/g2d_mirror/Streptococcus_pneumoniae_R6/ASM704v1_genomic.fna	
		/data/bagel4/bagel4_AOI_identification_hmm_rules.pl -dna /data/bagel4/examples/bottromycin.fna
" ;

&parseparam() ;

# ---------------------------------------------------------- main -------------------------------------------------------------------------------


# testing:
# cd /tmp/BAGEL4WRAPPER/129.125.143.35gq19j8k9b8f1rsjq94v6bml6o5492
# sudo -u www-data /data/bagel4/bagel4_AOI_identification_hmm_rules_v2.pl -s /tmp/BAGEL4WRAPPER/129.125.143.35gq19j8k9b8f1rsjq94v6bml6o5492 -dna /tmp/BAGEL4WRAPPER/129.125.143.35gq19j8k9b8f1rsjq94v6bml6o5492/AP0072091.0.fna -queryname AP0072091.0


# 1. translate in all 6 frames
	my %sixframes = anne_genomics::fasta2hash("$sessiondir/$queryname.six_frames.fna") ;
	my $frame_length = length($sixframes{frame_1}) ;
	my $DNA_length = 3*$frame_length ;

	
# 2. Load the AOI identification rules
	my %PFAM_rules = get_all_hmm_rules() ;  # number, name, array of hmms
	my @HMMs = get_all_pfams();


# 3. Here the search starts for proximity HMMs:
	my %HMM_results_table = search_HMMs("$sessiondir/$queryname.six_frames.fna" , @HMMs) ;
	six_frames_PFAM_to_DNA_pos() ;


# 4. Check the HMMrules
	my %valid_rules = apply_rules() ;

	
# 5. Export the result table
	report_v2();

	
	
# -----------------------------------------------------------  functions  ------------------------------------------------------------------------
 
	sub show_HMM_results_table {
		# dummy subroutine for debugging purposes
		print "HMM_results_table\n" ;
		foreach my $ID (keys %HMM_results_table) { 
			print "$HMM_results_table{$ID}{Name}\t";	
			print "$HMM_results_table{$ID}{PFAM_name}\t";	
			print "$HMM_results_table{$ID}{PFAM_ID}\t";
			print "$HMM_results_table{$ID}{Evalue}\t";	
			print "$HMM_results_table{$ID}{Score}\t";	
			print "$HMM_results_table{$ID}{start}\t";	
			print "$HMM_results_table{$ID}{end}\n";
		}	
	}	


	sub get_all_hmm_rules {
		# Read the rules table and return the table in a hash
		print "Reading HMM rules\n" ;
		my @lines = anne_files::read_lines($AOI_identification_rules_table) ;
		my %result ;
		my $i = 0 ;
		foreach my $line (@lines) {
			my @items = split "\t", $line ;
			if (scalar @items >1) {
				$i++ ;
				$result{$i}{line} = $line ; 
				$result{$i}{Name} = shift @items ; 
				@{$result{$i}{hmms}} = @items ; 
				# print "\t$result{$i}{Name}\t";
				# print join ',',@{$result{$i}{hmms}} ;
				# print "\n";
			}
		}
		return %result ;
	}
	
 	sub get_all_pfams {
		# return all the HMMs of the Rules Table
		my @result ;
		foreach my $i (keys %PFAM_rules) {
			push (@result, @{$PFAM_rules{$i}{hmms}}) ; 
		}
		return anne_files::unique_array(@result) ;
	}

	
	
	sub hmm_inrange {
		my $key = shift ;  # the key of the rules
		my $FirstHMM = $PFAM_rules{$key}{hmms}[0] ;
		my %positions ;
		
		foreach my $ID (keys %HMM_results_table) {  # get all the positions of the first hmm of the rules
			if ( $HMM_results_table{$ID}{PFAM_ID} =~ m/$FirstHMM/ ) { 
				my $pos = $HMM_results_table{$ID}{dnastart}; 
				$positions{$pos} = 1 ;
			}
		} 
		
		for (my $i=1; $i<=(scalar(@{$PFAM_rules{$key}{hmms}}))-1; $i++) {
			my $hmm = $PFAM_rules{$key}{hmms}[$i] ;
			#print "\t\tCheck in_range $hmm to $FirstHMM\n";
			foreach my $pos (sort keys %positions) {
				foreach my $ID (keys %HMM_results_table) { 
					if ( $HMM_results_table{$ID}{PFAM_ID} =~ m/$hmm/ ) {
						my $dist = abs($pos-$HMM_results_table{$ID}{dnastart}) ; 
						if ($dist > 0 and $dist < $distance ) { 
							print "\t\t\t Range: pos_primary=$pos hmm=$FirstHMM poshmm=$HMM_results_table{$ID}{dnastart} distance=$dist\n" ;
							$positions{$pos}++ ; # add a hit of the hmm is found in RAnge
							last ;
						}
					}	
				}
			}	
		}
		my @result  ;
		foreach my $pos (sort keys %positions) {  # return primary positions if in range with all HMMs of the rule
			if ($positions{$pos} == scalar(@{$PFAM_rules{$key}{hmms}}) ) { push @result, $pos; } 
		}	
		return @result ;
	}
	
	
	sub apply_rules {
		my %result ;  # contain the positions of valid rules
		my $resultkey = 0 ;
		print "Checking Rules:\n" ; 
		foreach my $key (sort keys %PFAM_rules) {
			# check if rules are valid e.g. PF05147	PF04055
			print "\t$PFAM_rules{$key}{line}\n";
			#print "\tNumber of hmms ==>".(scalar(@{$PFAM_rules{$key}{hmms}}))."\n";
			if (scalar @{$PFAM_rules{$key}{hmms}} > 1) { 
				my $hmm = $PFAM_rules{$key}{hmms}[0] ;
				my @positions ;
				my @valid_positions = hmm_inrange($key) ;
				foreach my $valid_position (@valid_positions) { 
					print  "\t\tMulti rule ==> $PFAM_rules{$key}{line} $hmm TRUE\n"; 
					$resultkey++ ;
					$result{$resultkey}{PFAM_rules} = $key ;
					$result{$resultkey}{Position} = $valid_position ;
				}	

			} else {
				# only one HMM, just get the result
				my $hmm = $PFAM_rules{$key}{hmms}[0] ;
				#print "\tSingle rule ==>$PFAM_rules{$key}{line}\n";
				foreach my $ID (keys %HMM_results_table) { 
					if ( $HMM_results_table{$ID}{PFAM_ID} =~ m/$hmm/ ) { 
						$resultkey++ ;
						print "\t$hmm\t".$HMM_results_table{$ID}{dnastart}."\tTRUE\n"; 
						$result{$resultkey}{PFAM_rules} = $key ;
						$result{$resultkey}{Position} = $HMM_results_table{$ID}{dnastart} ;
					}
				}
			}	
		}
		return %result ;
	} 
	
 
 
sub report_v2 {
	print "\n\n------------------ report bagel4_AOI_identification_hmm_rules --------------------------\n";
	print "DNA file 	= $fna_file\n";
	print "DNA length 	= $DNA_length\n";
	my @table ;
	my @fna ;
	my $DNA = anne_genomics::get_dna_from_fasta($fna_file) ;
	my $count = 0 ;
	foreach my $key (sort keys %valid_rules) {
		$count++ ;
		my $RecordName= "$queryname.primaryHMM.$count" ;
		my $start = $valid_rules{$key}{Position}-($distance/2) ;
		my $end = $start + $distance ;
		my $RuleKey = $valid_rules{$key}{PFAM_rules} ;
		push @fna, ">$RecordName" ;
		push @fna, substr $DNA, $start , $distance ;	
		push @table, "$RecordName\t$start\t$end\t+\t$PFAM_rules{$RuleKey}{Name}\thmm";
	}
	anne_files::write_lines("$sessiondir/$queryname.AOI.hmm.table", @table);
	anne_files::write_lines("$sessiondir/$queryname.AOI.hmm.fna", @fna);
	print "\n";
	print "AOIs found on the basis of HMM rules:\n------------------------------------------------------------\n".(join "\n", @table)."\n" ;
	print "Result Table written to $queryname.AOI.hmm.table\n" ;
} 


sub six_frames_PFAM_to_DNA_pos {
	# convert the protein position to the DNA position
	#my $DNA_length = length(anne_genomics::get_dna_from_fasta($fna_file)) ;
	foreach my $ID (keys %HMM_results_table) {
		foreach my $item (keys %{$HMM_results_table{$ID}}) {
			if ($HMM_results_table{$ID}{Name} =~ m/^frame_(\d)/) {
				$HMM_results_table{$ID}{dnastart} = $HMM_results_table{$ID}{start} * 3  ;
				$HMM_results_table{$ID}{dnaend}   = $HMM_results_table{$ID}{end} * 3 ;
				$HMM_results_table{$ID}{strand} = '+' ;
				if ($1=='4' or $1=='5' or $1=='6') {   # frame 4, 5 and 6 are from the lower strand
					$HMM_results_table{$ID}{dnastart} = $DNA_length - $HMM_results_table{$ID}{dnastart} ;
					$HMM_results_table{$ID}{dnaend}   = $DNA_length - $HMM_results_table{$ID}{dnaend} ;
					$HMM_results_table{$ID}{strand} = '-' ;
				}
			}
		}	
	}
}



sub hmmsearch {
	my ($protein_fasta_file, @HMMs) = @_ ;
	my @result_files ;
	my @log ;
	my $domtbloutDir = "$sessiondir/hmm_rules_domtblout" ;
	mkdir $domtbloutDir, 0755 if ( !-d $domtbloutDir );
	foreach my $pfam ( @HMMs ) {
		print "Screening for $pfam\n" ;
		my $domtblout = "$domtbloutDir/$queryname.$pfam.domtblout";
		my $tmp = "$conf{hmmsearch}/hmmsearch --noali --cpu $cpu --domtblout $domtblout --domT $domT $program_dir/$conf{HMMs_folder	}/$pfam.hmm $protein_fasta_file >>$domtbloutDir/logfile.txt" ;
		push @log, $tmp ;
		system($tmp) ;
		#push (@result_files, $domtblout) if (check_hmm_found($domtblout));
		push (@result_files, $domtblout);
	}
	anne_files::write_lines("$domtbloutDir/00.hmmsearch.log", @log) ;
	return @result_files ;
}

sub search_HMMs {
	my ($protein_fasta_file, @HMMs) = @_ ;	# fasta file of proteins and the HMMs to be searched
	my @result_files = hmmsearch($protein_fasta_file, @HMMs) ; 

	# parse the results to a table
	my $ID = 0 ;
	my %result_table ;
	#print "\nHMMs found:\n";
	foreach my $result_file (@result_files) {
		my @lines = anne_files::read_lines($result_file);
		foreach my $line (@lines) {
			if ($line =~ m/^frame_/) {
				$ID++ ;
				my @items = split /\s+/, $line ;
				$result_table{$ID}{Name} 		= $items[0] ;
				$result_table{$ID}{PFAM_name}	= $items[3] ;
				$result_table{$ID}{PFAM_ID} 	= $items[4] ;
				$result_table{$ID}{Evalue} 		= $items[6] ;
				$result_table{$ID}{Score} 		= $items[13] ;
				$result_table{$ID}{start} 		= $items[17] ;
				$result_table{$ID}{end} 		= $items[18] ;
				#print "\t".$items[0]."\t".$items[3]."\t".$items[4]."\t".$items[6]."\t".$items[13]."\t".$items[17]."\t".$items[18]."\n" ;
			}
		}	
	}

	return %result_table ;
}	
	
sub parseparam {
    my $var ;
    my @arg = @ARGV ;
    while(@arg) {
        $var = shift(@arg) ;
        die $usage if ($var eq '-h' or $var eq '--help') ;
		$sessiondir		= shift(@arg) if($var eq '-s') ;
		$fna_file		= shift(@arg) if($var eq '-dna') ;
		$queryname		= shift(@arg) if($var eq '-queryname') ;
		$cpu			= shift(@arg) if($var eq '-cpu') ;
		$contextsize	= shift(@arg) if($var eq '-cz') ;
		$outputfile		= shift(@arg) if($var eq '-o') ;
		$debug		= shift(@arg) if($var eq '-debug') ;
    }
    die $usage if (!$fna_file) ;
}


