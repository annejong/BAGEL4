#!/usr/bin/env perl

use strict ;
use warnings ;
use lib "/data/bagel4/lib" ;
use bagel4_functions ;
use lib "/data/molgentools/lib";
use anne_files ;
use File::Basename;


# ---------------------------------------------------------- parameters -------------------------------------------------------------------------

my $sessiondir = '/usr/bagel4/test';
my $program_dir = dirname($0) ;
my $AOI_identification_rules_table = "$program_dir/tables/AOI_identification_rules.txt" ;
my $outputfile = 'hmmsearch.results';
my $queryname = 'queryname';
my $fna_file ;
my $debug = 1 ;

my %conf = bagel4_functions::read_conf("$program_dir/bagel4.conf") ;
my $contextsize = $conf{contextsize};
my $domT = $conf{domT} ; # PFAM score cutoff --domT
my $cpu = $conf{cpu} ;

my $usage = "options:
	-s Sessiondir [default=$sessiondir]
	-dna dna input file in FASTA format
	-cpu number of cpu's used [default=$cpu]
	-queryname name of the query
	-size contextsize [default=$contextsize]
	-o Results file [default=$outputfile]
	
	this routine uses the identification rules from table AOI_identification_rules.txt [default=$AOI_identification_rules_table]
				
		e.g.  /usr/bagel4/bagel4_AOI_identification_hmm_rules.pl -dna /var/genomes/g2d_mirror/Streptococcus_pneumoniae_R6/ASM704v1_genomic.fna	
		/usr/bagel4/bagel4_AOI_identification_hmm_rules.pl -dna /usr/bagel4/examples/bottromycin.fna
" ;

&parseparam() ;

# ---------------------------------------------------------- main -------------------------------------------------------------------------------


# 1. translate in all 6 frames
	my %sixframes = anne_genomics::fasta2hash("$sessiondir/$queryname.six_frames.fna") ;
	my $frame_length = length($sixframes{frame_1}) ;
	my $DNA_length = 3*$frame_length ;

# 2. Load the AOI identification rules
	my %PFAM_rules = anne_files::read_table_to_hash($AOI_identification_rules_table) ;
	foreach my $rule_key (keys %PFAM_rules) { $PFAM_rules{$rule_key}{rule} =~ s/\(|\)//g ; } # clean the table by removing brackets

# 3. search all HMM profiles from the table tables/context_hmm.txt
	my @AOI_identification_HMMs ;
	foreach my $key ( keys %PFAM_rules ) {	push @AOI_identification_HMMs, $PFAM_rules{$key}{proximity} ; } # get the HMMs from the AOI_identification_rules_table
	@AOI_identification_HMMs = anne_files::unique_array(@AOI_identification_HMMs) ; # remove replicates

	foreach my $key (@AOI_identification_HMMs) { print "===> AOI_identification_HMM proximity: $key\n"; }

	
	# Here the search starts for proximity HMMs:
	my %HMM_results_table = search_HMMs("$sessiondir/$queryname.six_frames.fna" , @AOI_identification_HMMs) ;
	#write_PFAM_results("$sessiondir/$queryname.$outputfile");
	
# 4. get the protein AOI regions from the HMM positions
	my %AOI = protein_AOI_from_PFAM_positions() ; 

	
	
	
# 5. AOI indentification rules
	my %valid_AOI = new_apply_rules() ;
	#if ($debug) {%valid_AOI = new_apply_rules(); } else {%valid_AOI = apply_rules();}	
	
# 6. Export the complete result table
	report();

# -----------------------------------------------------------  functions  ------------------------------------------------------------------------
 
 
sub report {
	print "\n\n------------------ report bagel4_AOI_identification_hmm_rules --------------------------\n";
	print "DNA file 	= $fna_file\n";
	print "DNA length 	= $DNA_length\n";
	my @table ;
	my @fna ;
	my $DNA = anne_genomics::get_dna_from_fasta($fna_file) ;
	foreach my $key (sort keys %valid_AOI) {
		my $AOI_key  = $valid_AOI{$key}{AOI} ;
		my $rule_key = $valid_AOI{$key}{rule} ;
		if ($AOI{$AOI_key}{Name} =~ m/^frame_(\d)/) {
			$AOI{$AOI_key}{dnastart} = $AOI{$AOI_key}{protstart} * 3  ;
			$AOI{$AOI_key}{dnaend}   = $AOI{$AOI_key}{protend} * 3 ;
			$AOI{$AOI_key}{strand}   = '+';
			if ($1>3) {  # frame 4,5,6 from lower strand, recalculate from end
				$AOI{$AOI_key}{dnastart} = $DNA_length - $AOI{$AOI_key}{dnastart} ;
				$AOI{$AOI_key}{dnaend}   = $DNA_length - $AOI{$AOI_key}{dnaend} ;  
				$AOI{$AOI_key}{strand}   = '-';
			}
		
		}	
		if ($AOI{$AOI_key}{dnastart} > $AOI{$AOI_key}{dnaend}) {
			my $tmp = $AOI{$AOI_key}{dnastart} ; $AOI{$AOI_key}{dnastart} = $AOI{$AOI_key}{dnaend}; $AOI{$AOI_key}{dnaend} = $tmp ;
		}	
		print "\tPoximity PFAM = $AOI{$AOI_key}{centralPFAM} in $AOI{$AOI_key}{Name}\n";
		print "\tAOI_$AOI_key identified on the basis of rule ($PFAM_rules{$rule_key}{rule}) for $PFAM_rules{$rule_key}{name}\n";
		print "\tAOI frame from $AOI{$AOI_key}{protstart} to $AOI{$AOI_key}{protend}\n";
		print "\tAOI DNA   from $AOI{$AOI_key}{dnastart} to $AOI{$AOI_key}{dnaend}\n";
		push @table, "$queryname.AOI_$AOI_key\t$AOI{$AOI_key}{dnastart}\t$AOI{$AOI_key}{dnaend}\t$AOI{$AOI_key}{strand}\t$PFAM_rules{$rule_key}{name}\thmm";
		my $len = abs($AOI{$AOI_key}{dnaend} - $AOI{$AOI_key}{dnastart});
		print "==============> Length=$len\n";
		print "==============> Start=$AOI{$AOI_key}{dnastart}\n";
		push @fna, ">$queryname.AOI_$AOI_key" ;
		push @fna, substr $DNA, $AOI{$AOI_key}{dnastart}, $len ;
	}
	anne_files::write_lines("$sessiondir/$queryname.AOI.hmm.table", @table);
	anne_files::write_lines("$sessiondir/$queryname.AOI.hmm.fna", @fna);
	print "AOI.hmm.table written to $queryname.AOI.hmm.table\n" ;
}


sub get_rules_from_centralPFAM {
	my $pfam = shift ;
	my @results ;
	foreach my $key (keys %PFAM_rules) {
		if ($PFAM_rules{$key}{rule} =~ m/$pfam/) { push @results, $key ; }	
	}
	return @results ;
}

sub get_pfams_from_rules {
	my @results ;
	foreach my $key (@_) {
		my $rule = $PFAM_rules{$key}{rule} ;
		$rule =~ s/\(|\)//g ;  # remove brackets
		$rule =~ s/ AND /|/g ; # make separator
		my @items = split /\|/, $rule ;
		shift @items ;  
		push @results, @items ; # add the pfams except the first, because this is the proximity pfam
	}
	#print "PFAMs from rules:".(join ";", @results)."\n";
	return anne_files::unique_array(@results) ;
}


sub apply_rules {
	# here we check the rules for the initial AOIs
	my $AOIfilename = "$sessiondir/$queryname.tmp_proteins.fna" ;
	my %results ;
	my $result_key = 0 ;
	foreach my $key (keys %AOI) {
		print "\nChecking AOI $key for PFAM rules\n";
		# write all 6 frames of AOI to file. The secondary PFAM can also be in one of the other frames
		my @proteins ;
		for (my $i=1; $i<=6; $i++) { 
			my $frame = $i ;
			push @proteins, ">AOIprotein$frame";;
			push @proteins, substr $sixframes{'frame_'.$i},   $AOI{$key}{protstart}, ($AOI{$key}{protend}- $AOI{$key}{protstart}) ;
		}
		print "\tWrite 6 frame translation to $AOIfilename\n";
		anne_files::write_lines($AOIfilename, @proteins) ;
		
		# check the presence of the secondary PFAMs
		print "\tGet rules for AOI $AOI{$key}{Name} containing $AOI{$key}{centralPFAM}\n";
		
		my @rule_keys = get_rules_from_centralPFAM($AOI{$key}{centralPFAM}); # these rules should be checked
		my @secondpfams = get_pfams_from_rules(@rule_keys) ;
		print "\tSecondairy PFAMs from rules: ".(join ";", @secondpfams)."\n";
		my @secondpfams_found ;
		foreach my $secondpfam (@secondpfams) { 
			my @result_files = hmmsearch($AOIfilename, $secondpfam) ; # Only send one pfam to screen, so only one result file will be returned
			# /// remove from $secondpfam without hit
			my @lines = anne_files::read_lines($result_files[0]);
			# print "$lines[3]\n";
			if ($lines[3] =~ m/^AOIprotein/) { 
				print "\t\t\t>$secondpfam found in AOI $key\n";
				push @secondpfams_found, $secondpfam ; 
			} 
		}	
		
		# now check the rules
		foreach my $rule_key (@rule_keys) {
			# print "\tApplying rule $PFAM_rules{$rule_key}{rule} for $PFAM_rules{$rule_key}{name} \n";
			$PFAM_rules{$rule_key}{rule} =~ s/\(|\)//g ;
			my @ANDrules = split " AND ", $PFAM_rules{$rule_key}{rule} ;
			shift @ANDrules ;  # remove the first rule which is the proximity PFAM
			my $true_count = 0 ;
			foreach my $ANDrule (@ANDrules) {
				my @ORrules = split /\|/, $ANDrule ;
				foreach my $ORrule (@ORrules) { 
					if  (grep { /$ORrule/ } @secondpfams_found) { $true_count++ ; last; }
				}
			}
			if ($true_count == scalar @ANDrules) {
				print "\t===> AOI $key valid on the basis of rule $PFAM_rules{$rule_key}{rule} <===\n";
				$result_key++;
				$results{$result_key}{AOI} = $key ;
				$results{$result_key}{rule} = $rule_key ;
			}	
		}
	}
	return %results ;

}

sub system_command {
	my $command = shift ;
	system($command) ;
	anne_files::append_lines("$sessiondir/system_command_AOI_identification_hmm_rules.log", $command) ;
}

sub new_apply_rules {
	# here we check the rules for the initial AOIs
	my $AOIfilename = "$sessiondir/$queryname.tmp_proteins.fna" ;
	my %results ;
	my $result_key = 0 ;
	# check each AOI for the RULES
	foreach my $key (keys %AOI) {
		print "\nChecking AOI $key for HMM rules\n";
		# 1. Write all 6 frames of AOI to file. The secondary HMM can also be in one of the other frames
			my @proteins ;
			for (my $i=1; $i<=6; $i++) { 
				my $frame = $i ;
				push @proteins, ">AOIprotein$frame";;
				push @proteins, substr $sixframes{'frame_'.$i},   $AOI{$key}{protstart}, ($AOI{$key}{protend}- $AOI{$key}{protstart}) ;
			}
			print "\tWrite 6 frame translation to $AOIfilename\n";
			anne_files::write_lines($AOIfilename, @proteins) ;
		
		# 2. Check the presence of the secondary PFAMs for this Proximity HMM
			print "\tGet rules for AOI $AOI{$key}{Name} proximity $AOI{$key}{centralPFAM}\n";
			my @rule_keys = get_rules_from_centralPFAM($AOI{$key}{centralPFAM}); # these/this (multiple) rule(s) should be checked
			foreach my $rule_key (@rule_keys) {
				print "\t\tChecking rule_key $rule_key rule:$PFAM_rules{$rule_key}{rule}\n";
				my $rule = $PFAM_rules{$rule_key}{rule} ;
				$rule =~ s/ AND /|/g ; # make separator
				my @rule_pfams = split /\|/, $rule ;
				shift @rule_pfams ;  # add the pfams except the first, because this is the proximity pfam
				print "\t\tSecondairy PFAMs from rules: ".(join ";", @rule_pfams)."\n";	

				# do the hmmsearch here: instead of result files
				my @pfams_found ;				
				foreach my $pfam (@rule_pfams) {
					my $domtblout = "$sessiondir/$queryname.$pfam.domtblout";
					system_command("$conf{hmmsearch}/hmmsearch --noali --cpu $cpu --domtblout $domtblout --domT $domT $program_dir/$conf{HMMs_folder}/$pfam.hmm $AOIfilename >>$sessiondir/logfile.txt") ;
					my @lines = anne_files::read_lines($domtblout);
					if ($lines[3] =~ m/^AOIprotein/) { 
						print "\t\t-->$pfam found in AOI $key\n";
						push @pfams_found, $pfam ; 
					} 
				}

				# now check the rule
				# print "\tApplying rule $PFAM_rules{$rule_key}{rule} for $PFAM_rules{$rule_key}{name} \n";
				my @ANDrules = split " AND ", $PFAM_rules{$rule_key}{rule} ;
				shift @ANDrules ;  # remove the first rule which is the proximity PFAM
				my $true_count = 0 ;
				foreach my $ANDrule (@ANDrules) {
					my @ORrules = split /\|/, $ANDrule ;
					foreach my $ORrule (@ORrules) { 
						if  (grep { /$ORrule/ } @pfams_found) { $true_count++ ; last; }
					}
				}
				if ($true_count == scalar @ANDrules) {
					print "\t===> AOI $key valid on the basis of rule $PFAM_rules{$rule_key}{name}: $PFAM_rules{$rule_key}{rule} <===\n";
					$result_key++;
					$results{$result_key}{AOI} = $key ;
					$results{$result_key}{rule} = $rule_key ;
				}
			}
			
		}
	return %results ;
}

sub protein_AOI_from_PFAM_positions {
	my %AOI ;
	my $key = 0 ;
	foreach my $ID (keys %HMM_results_table) {
		$HMM_results_table{$ID}{PFAM_ID} =~ s/\..*// ;  # remove version number
		$key++ ;
		$AOI{$key}{centralPFAM} = $HMM_results_table{$ID}{PFAM_ID} ;
		$AOI{$key}{strand} 		= $HMM_results_table{$ID}{strand} ;
		$AOI{$key}{Name} 		= $HMM_results_table{$ID}{Name} ;
		$AOI{$key}{frame} 		= $HMM_results_table{$ID}{Name} ;
		$AOI{$key}{frame}  		=~ s/frame_//;
		$AOI{$key}{protstart} 	= $HMM_results_table{$ID}{start} - int($contextsize/6) ;
		$AOI{$key}{protstart} 	= 1 if ($AOI{$key}{protstart}<1) ;
		$AOI{$key}{protend} 	= $HMM_results_table{$ID}{end} + int($contextsize/6) ;
		$AOI{$key}{protend} 	= $frame_length if ($AOI{$key}{protend}>$frame_length) ;
	}	
	return %AOI ;
}

sub get_PFAMs {
	my $rule = shift ;
	$rule =~ s/\(|\)//g ;
	return split " AND ", $rule ;
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
	my @debug ;
	foreach my $pfam ( @HMMs ) {
		my $domtblout = "$sessiondir/$queryname.$pfam.domtblout";
		my $tmp = "$conf{hmmsearch}/hmmsearch --noali --cpu $cpu --domtblout $domtblout --domT $domT $program_dir/$conf{HMMs_folder	}/$pfam.hmm $protein_fasta_file >>$sessiondir/logfile.txt" ;
		push @debug, $tmp ;
		system($tmp) ;
		#push (@result_files, $domtblout) if (check_hmm_found($domtblout));
		push (@result_files, $domtblout);
	}
	anne_files::write_lines("$sessiondir/00.debug", @debug) ;
	return @result_files ;
}

sub search_HMMs {
	my ($protein_fasta_file, @HMMs) = @_ ;	# fasta file of proteins and the HMMs to be searched
	my @result_files = hmmsearch($protein_fasta_file, @HMMs) ; 

	# parse the results to a table
	my $ID = 0 ;
	my %result_table ;
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
				print "HMM found:\t".$items[0]."\t".$items[3]."\t".$items[4]."\t".$items[6]."\t".$items[13]."\t".$items[17]."\t".$items[18]."\n" ;
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


