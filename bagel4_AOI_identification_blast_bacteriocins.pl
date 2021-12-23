#!/usr/bin/env perl

use strict ;
use warnings ;
use lib "/data/bagel4/lib" ;
use bagel4_functions ;
use lib "/data/molgentools/lib";
use anne_files ;
use anne_genomics ;
use File::Basename;


# ---------------------------------------------------------- parameters -------------------------------------------------------------------------

my $sessiondir = '.';
my $program_dir = dirname($0) ;
my $translate = 'false';
my $queryname = 'query';
my $fna_file = '/data/g2d_mirror/Streptococcus_pneumoniae_R6/ASM704v1_genomic.fna' ;
my $outputfile = "AOI.blast.fna";

my %conf = bagel4_functions::read_conf("$program_dir/bagel4.conf") ;
my $contextsize 	= $conf{contextsize};
my $cpu = $conf{cpu} ;

my $usage = "option:
		-s		sessiondir [default=$sessiondir]
		-dna	dna input file in FASTA format
		-size	contextsize [default=$contextsize]
		-queryname name of the query
		-translate false|true [default=$translate], translation is already done for the hmm search
		-o		results file [default=$outputfile]
				
		this routine will blast all six frames to the 3 bacteriocin databases for AOI identification
				
		e.g.  time /usr/bagel4/bagel4_blast_bacteriocins.pl -translate true -dna /var/genomes/g2d_mirror/Streptococcus_pneumoniae_R6/ASM704v1_genomic.fna	
		/usr/bagel4/bagel4_blast_context_search.pl -translate true -dna /usr/bagel4/examples/bottromycin.fna
" ;

&parseparam() ;
my $logfile = "$sessiondir/log_blast_AOI.txt";

# ---------------------------------------------------------- main -------------------------------------------------------------------------------

# 1. Container for results
	my %AOI ;
	my $AOI_key = 0 ;
	my $DNA = anne_genomics::get_dna_from_fasta($fna_file) ;
	my $DNA_length = length($DNA) ;
	my %six_frames = anne_genomics::fasta2hash("$sessiondir/$queryname.six_frames.fna");
	
# 2. Blast all ORFs against bacteriocin database I, II and III on all six frames
	anne_files::write_log($logfile,"\n------------ step 2 Blast six frames against Bacteriocin databases------------\n",'true');
	my $command = "formatdb -p T -o T -i $sessiondir/$queryname.six_frames.fna -l $sessiondir/formatdb.log"; system($command) ;
	my %blast = blast_frames("$program_dir/$conf{bacteriocinI_proteins}", $conf{blast_evalue_bacteriocinI}, 'blast_bacteriocinI.results');
	my @hits ;
	if (scalar (keys %blast)) { push @hits, AOI_identification('bacteriocinI'); }
	%blast = blast_frames("$program_dir/$conf{bacteriocinII_proteins}", $conf{blast_evalue_bacteriocinII},  'blast_bacteriocinII.results');
	if (scalar (keys %blast)) { push @hits, AOI_identification('bacteriocinII'); }
	%blast = blast_frames("$program_dir/$conf{bacteriocinIII_proteins}", $conf{blast_evalue_bacteriocinIII},'blast_bacteriocinIII.results');
	if (scalar (keys %blast)) { push @hits, AOI_identification('bacteriocinIII'); }
	
# 3. Export AOI fasta file
	anne_files::write_lines("$sessiondir/$queryname.AOI.blast.table", @hits);
	
	my @lines ;
	foreach my $AOI_key (keys %AOI) {
		push @lines, ">$queryname.AOI_$AOI_key";
		my $len = abs($AOI{$AOI_key}{dnaend}-$AOI{$AOI_key}{dnastart});
		push @lines, substr $DNA, $AOI{$AOI_key}{dnastart}, $len ;
		#print "writing >$queryname.AOI_$AOI_key $AOI{$AOI_key}{dnastart}, $AOI{$AOI_key}{dnaend} $AOI{$AOI_key}{strand}\t$len bp\n";
	}
	anne_files::write_lines("$sessiondir/$queryname.$outputfile", @lines) ;
	print "AOI DNA sequences written to $sessiondir/$queryname.$outputfile\n";
	print "AOI Blast table   written to $sessiondir/$queryname.AOI.blast.table\n";
	
# -----------------------------------------------------------  functions  ------------------------------------------------------------------------


sub AOI_identification {
	my $class = shift ;
	print "AOI identification on the basis of blast hit found in $class database\n"; 
	my @lines ;
	foreach my $key (sort keys %blast) {
		$AOI_key++ ;
		$AOI{$AOI_key}{class} = $class ;
		if ($blast{$key}{Subject} =~ m/^frame_(\d)/) {
			$AOI{$AOI_key}{dnastart} = $blast{$key}{s_start} * 3 - int($contextsize/2) ;
			$AOI{$AOI_key}{dnastart} = 1 if ($AOI{$AOI_key}{dnastart} < 1) ;
			$AOI{$AOI_key}{dnaend}   = $blast{$key}{s_end} * 3   + int($contextsize/2) ;
			$AOI{$AOI_key}{dnaend}   = $DNA_length if ($AOI{$AOI_key}{dnaend}>$DNA_length) ;
			$AOI{$AOI_key}{strand}   = '+'; 
			$AOI{$AOI_key}{bit_score}   = $blast{$key}{bit_score}; 
			if ($1>3) {  # frame 4,5,6 from lower strand, recalculate from end
				my $tmp_start = $AOI{$AOI_key}{dnastart}  ;
				$AOI{$AOI_key}{dnastart}   = $DNA_length - $AOI{$AOI_key}{dnaend} ;
				$AOI{$AOI_key}{dnaend} = $DNA_length - $tmp_start ;
				$AOI{$AOI_key}{strand}   = '-'; 
			}
			
		}
		print "\t$queryname.AOI_$AOI_key ($AOI{$AOI_key}{dnastart} - $AOI{$AOI_key}{dnaend}) identified on the basis of blast hit with $blast{$key}{Query}\n";
		if ($AOI{$AOI_key}{strand} eq '+') {
			push @lines, "$queryname.AOI_$AOI_key\t$AOI{$AOI_key}{dnastart}\t$AOI{$AOI_key}{dnaend}\t$AOI{$AOI_key}{strand}\t$blast{$key}{Query}\t$class\t$AOI{$AOI_key}{bit_score}"; 
		} else {
			push @lines, "$queryname.AOI_$AOI_key\t$AOI{$AOI_key}{dnaend}\t$AOI{$AOI_key}{dnastart}\t$AOI{$AOI_key}{strand}\t$blast{$key}{Query}\t$class\t$AOI{$AOI_key}{bit_score}"; 
		}	
	}
	return @lines ;
	
}

sub blast_frames {
	my ($db, $evalue, $blastresultfile) = @_ ;  # $resultlable is the label in the result_table has file
	$blastresultfile = "$sessiondir/$queryname.$blastresultfile";
	# swapped -db and -query to speed up the blast session
	my $tmp = "blastp -outfmt 6 -query $db -db $sessiondir/$queryname.six_frames.fna -max_target_seqs 1 -num_threads $cpu -evalue $evalue -out $blastresultfile" ;
	system($tmp) ;
	&anne_files::add_header($blastresultfile, "Query\tSubject\tpercent\talign_len\tmistmatches\tgap\tq_start\tq_end\ts_start\ts_end\tevalue\tbit_score");
	my %blast = anne_files::read_table_to_hash($blastresultfile) ;
	
	# check here if the cross blast also gives a hit and remove the hit if not
	# Query	Subject	percent	align_len	mistmatches	gap	q_start	q_end	s_start	s_end	evalue	bit_score
	# 119.2;Lactococcin_972_(Lcn972)	frame_3	34.34	99	55	4	1	90	206588	206685	1e-08	50.8
	foreach my $key (keys %blast) {
		my $frame = $blast{$key}{Subject} ;
		my $protein = substr $six_frames{$frame}, $blast{$key}{s_start}, $blast{$key}{s_end}-$blast{$key}{s_start} ;
	 	anne_files::write_string("$sessiondir/tmp.faa", ">$key\n$protein\n") ;
		my $command = "blastp -outfmt 6 -db $db -query $sessiondir/tmp.faa -max_target_seqs 1 -num_threads $cpu -evalue $evalue -out $sessiondir/tmp.blast" ;
		#print "===>$command\n";
		system($command) ;
	 	if (-z "$sessiondir/tmp.blast") { delete($blast{$key}) ; print "===========> deleting $key <========================\n" ;} # not hit found if file is empty
	}
	return %blast ;
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
		$contextsize	= shift(@arg) if($var eq '-cz') ;
		$outputfile		= shift(@arg) if($var eq '-o') ;
    }
    die $usage if (!$fna_file) ;
}


