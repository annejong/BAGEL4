#!/usr/bin/env perl

use strict ;
use warnings ;
use bagel4_functions ;
use lib "/usr/molgentools/lib";
use anne_files ;
use File::Basename;


# ---------------------------------------------------------- parameters -------------------------------------------------------------------------

my $sessiondir = '/usr/bagel4/test';
my $program_dir = dirname($0) ;
my $contextsize = 20000 ;
my $outputfile = 'hmmsearch.results';
my $domT = 25 ; # PFAM score cutoff --domT
my $hmmsearch_cpu = 16 ;


my $fna_file = '/var/genomes/g2d_mirror/Streptococcus_pneumoniae_R6/ASM704v1_genomic.fna' ;
my $usage = "./bagel4_hmmsearch.pl
				-s Sessiondir [default=$sessiondir]
				-dna dna input file in FASTA format
				-cpu number of cpu's used [default=$hmmsearch_cpu]
				-size contextsize [default=$contextsize]
				-o Results file [default=$outputfile]
		e.g.  /usr/bagel4/bagel4_hmmsearch.pl -dna /var/genomes/g2d_mirror/Streptococcus_pneumoniae_R6/ASM704v1_genomic.fna	
" ;

&parseparam() ;

# ---------------------------------------------------------- main -------------------------------------------------------------------------------


# 1. translate in all 6 frames
	my %sixframes = translate_six_frames($fna_file) ;

# 2. Load the AOI identification rules
	my %PFAM_rules = anne_files::read_table_to_hash("$program_dir/tables/AOI_identification_rules.txt") ;

# 3. search all HMM profiles from the table tables/context_hmm.txt
	my @proximity_pfams ;
	foreach my $key ( keys %PFAM_rules ) {	push @proximity_pfams, $PFAM_rules{$key}{proximity} ; } # get the proximity_pfams from the AOI rules table
	@proximity_pfams = anne_files::unique_array(@proximity_pfams) ; # remove replicates
	my %PFAM_results_table = search_PFAMS_v4("$sessiondir/six_frames.fna" , @proximity_pfams) ;
	

	
# 4. calculate the DNA position and add this to %PFAM_results_table
	#six_frames_PFAM_to_DNA_pos() ;
		#		{dnastart}
		#		{dnaend}  
		#		{strand}
	my %AOI = protein_AOI_from_PFAM_positions() ; # get the protein AOI regions from the PFAM positions

# 5. AOI indentification rules
	apply_rules();	

# 6. Export the complete result table
#write_PFAM_results("$sessiondir/$outputfile");


# -----------------------------------------------------------  functions  ------------------------------------------------------------------------



sub get_rules_from_centralPFAM {
	my $pfam = shift ;
	my @results ;
	foreach my $key (keys %PFAM_rules) {
		if ($PFAM_rules{$key}{rule} =~ m/$pfam/) {
			print "\t\t".$PFAM_rules{$key}{name}."\t".$PFAM_rules{$key}{rule}."\n";
			push @results, $key ;
		}	
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
	my $AOIfilename = "$sessiondir/AOI.fna" ;
	foreach my $key (keys %AOI) {
		print "\nChecking AOI $key for PFAM rules\n";
		# write AOI to file
		my $prot = substr $sixframes{$AOI{$key}{frame}}, $AOI{$key}{protstart}, $AOI{$key}{protend} ;
		anne_files::write_string($AOIfilename, ">AOIprotein\n$prot\n") ;
		print "\tGet rules for AOI $AOI{$key}{Name} containing $AOI{$key}{centralPFAM}\n";
		my @rule_keys = get_rules_from_centralPFAM($AOI{$key}{centralPFAM});
		my @secondpfams = get_pfams_from_rules(@rule_keys) ;
		print "\tSecondairy PFAMs from rules: ".(join ";", @secondpfams)."\n";
		#my %secondPFAM_results = search_PFAMS_v4($AOIfilename, @secondpfams) ;
		foreach my $secondpfam (@secondpfams) { 
			my @result_files = hmmsearch($AOIfilename, $secondpfam) ; 
			# /// remove from $secondpfam without hit
			my @lines = anne_files::read_lines($result_files[0]);
			print "$lines[3]\n";
			if ($lines[3] =~ m/^AOIprotein/) { print "------------------------------------>Hit for $secondpfam\n"; } else  { print "NO hit for $secondpfam\n"; }
		}	
		
		# now check the rules
		foreach my $rule_key (@rule_keys) {
			print "\tApplying rule $PFAM_rules{$rule_key}{rule} for $PFAM_rules{$rule_key}{name} \n";
			$PFAM_rules{$rule_key}{rule} =~ s/\(|\)//g ;
			my @ANDrules = split " AND ", $PFAM_rules{$rule_key}{rule} ;
			shift @ANDrules ;  # remove the first rule which is the proximity PFAM
			my $true_count = 0 ;
			foreach my $ANDrule (@ANDrules) {
				my @ORrules = split /\|/, $ANDrule ;
				foreach my $ORrule (@ORrules) { 
					if  (grep { /$ORrule/ } @secondpfams) { $true_count++ ; last; }
				}
			}
			if ($true_count == scalar @ANDrules) {
				print "AOI valid!!\n";
			}	
		}
	}

}

sub protein_AOI_from_PFAM_positions {
	my %AOI ;
	my $key = 0 ;
	foreach my $ID (keys %PFAM_results_table) {
		$PFAM_results_table{$ID}{PFAM_ID} =~ s/\..*// ;  # remove version number
		$key++ ;
		$AOI{$key}{centralPFAM} = $PFAM_results_table{$ID}{PFAM_ID} ;
		$AOI{$key}{protstart} 	= $PFAM_results_table{$ID}{start} - int($contextsize/6) ;
		$AOI{$key}{protend} 	= $PFAM_results_table{$ID}{end} + int($contextsize/6) ;
		$AOI{$key}{strand} 		= $PFAM_results_table{$ID}{strand} ;
		$AOI{$key}{Name} 		= $PFAM_results_table{$ID}{Name} ;
		$AOI{$key}{frame} 		= $PFAM_results_table{$ID}{Name} ;
		$AOI{$key}{frame}  		=~ s/frame_//;
		#push @{$result{$PFAM_results_table{$ID}{PFAM_ID}}}, $PFAM_results_table{$ID}{dnastart} ; # this position is good enough . Better but slower would be between dnastart-dnaend
		print "protein AOI $key: on the basis of $PFAM_results_table{$ID}{PFAM_ID} in $AOI{$key}{Name} from $AOI{$key}{protstart} to $AOI{$key}{protend}\n";
	}	
	return %AOI ;
}

sub get_PFAMs {
	my $rule = shift ;
	$rule =~ s/\(|\)//g ;
	return split " AND ", $rule ;
}


sub write_PFAM_results {
	my $outputfilename = shift ;
	my @headers = ('Name','PFAM_name','PFAM_ID','Evalue','Score','start','end','dnastart','dnaend','strand') ; 	
	my @lines = join "\t", @headers ;
	foreach my $ID (sort { $PFAM_results_table{$a}{PFAM_name} cmp $PFAM_results_table{$b}{PFAM_name} || $PFAM_results_table{$a}{Name} cmp $PFAM_results_table{$b}{Name}} keys %PFAM_results_table) {
		my @row ;
		foreach my $header (@headers) {
			push @row, $PFAM_results_table{$ID}{$header}
		}	
		push @lines, join "\t", @row ;
	}
	#foreach my $line (@lines) { print "$line\n" ; }
	anne_files::write_lines($outputfilename, @lines) ;
}


sub six_frames_PFAM_to_DNA_pos {
	# convert the protein position to the DNA position
	my $DNA_len = length(anne_genomics::get_dna_from_fasta($fna_file)) ;
	foreach my $ID (keys %PFAM_results_table) {
		foreach my $item (keys $PFAM_results_table{$ID}) {
			if ($PFAM_results_table{$ID}{Name} =~ m/^frame_(\d)/) {
				$PFAM_results_table{$ID}{dnastart} = $PFAM_results_table{$ID}{start} * 3  ;
				$PFAM_results_table{$ID}{dnaend}   = $PFAM_results_table{$ID}{end} * 3 ;
				$PFAM_results_table{$ID}{strand} = '+' ;
				if ($1=='4' or $1=='5' or $1=='6') {   # frame 4, 5 and 6 are from the lower strand
					$PFAM_results_table{$ID}{dnastart} = $DNA_len - $PFAM_results_table{$ID}{dnastart} ;
					$PFAM_results_table{$ID}{dnaend}   = $DNA_len - $PFAM_results_table{$ID}{dnaend} ;
					$PFAM_results_table{$ID}{strand} = '-' ;
				}	
			}
		}	
	}
}


sub translate_six_frames {
	# translate DNA in all 6 frames and return the filename of the proteins
	my $FNA_file = shift;
	my @FNA = anne_files::read_lines($FNA_file) ;
	shift @FNA ; # remove header
	my $DNA = join('', @FNA) ;
	my %sixframes = bagel4_functions::six_frames($DNA) ;
	my @lines ;
	foreach my $frame (sort keys %sixframes) {
		push @lines, ">frame_$frame";
		push @lines, $sixframes{$frame} ;
	}	
	anne_files::write_lines("$sessiondir/six_frames.fna", @lines) ;
	return %sixframes ;
}

sub hmmsearch {
	my ($protein_fasta_file, @pfams) = @_ ;
	my @result_files ;
	foreach my $pfam ( @pfams ) {
		my $domtblout = "$sessiondir/$pfam.domtblout";
		push (@result_files, $domtblout) ;
		my $tmp = "/usr/bagel4/hmmsearch/hmmsearch/binaries/hmmsearch --noali --cpu $hmmsearch_cpu --domtblout $domtblout --domT $domT /usr/bagel4/hmm_db/$pfam.hmm $protein_fasta_file >>$sessiondir/logfile.txt" ;
		print "\t\tSearch $pfam in $protein_fasta_file\n";
		system($tmp) ;
	}
	return @result_files ;
}

sub search_PFAMS_v4 {
	my ($protein_fasta_file, @pfams) = @_ ;		# fasta file of proteins and the PFAMs to be searched
	my @result_files = hmmsearch($protein_fasta_file, @pfams) ; 

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
				#print $items[0]."\t".$items[3]."\t".$items[4]."\t".$items[6]."\t".$items[13]."\t".$items[17]."\t".$items[18]."\n" ;
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
		$hmmsearch_cpu	= shift(@arg) if($var eq '-cpu') ;
		$contextsize	= shift(@arg) if($var eq '-cz') ;
		$outputfile		= shift(@arg) if($var eq '-o') ;
    }
    die $usage if (!$fna_file) ;
}


