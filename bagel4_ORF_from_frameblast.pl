#!/usr/bin/env perl


# Uncover ORFs not found by small-orf or Glimmer in AOI_annotation




use strict ;
use warnings ;
use lib "/data/bagel4/lib" ;
use bagel4_functions ;
use lib "/data/molgentools/lib";
use anne_files ;
use anne_genomics ;
use File::Basename;


# ---------------------------------------------------------- parameters -------------------------------------------------------------------------

my $sessiondir = '/tmp/BAGEL4WRAPPER/test';
my $program_dir = dirname($0) ;
my $queryname  ;
my $usage = "option:
		-s		sessiondir [default=$sessiondir]
		-queryname	 

				
	e.g.  /data/bagel4/bagel4_ORF_from_frameblast.pl -s /tmp/BAGEL4WRAPPER/129.125.142.955g88tcaq2qoa0e0c73q73d6fu7863 -queryname CP0038401AOI_01.AOI_01	
	e.g.  /data/bagel4/bagel4_ORF_from_frameblast.pl -s /tmp/BAGEL4WRAPPER/129.125.142.955g88tcaq2qoa0e0c73q73d6fu7887 -queryname NC_0045672.AOI_01	
" ;
my %conf = bagel4_functions::read_conf("$program_dir/bagel4.conf") ;

&parseparam() ;

# ---------------------------------------------------------- main -------------------------------------------------------------------------

# 1. Translate to 6 frames
	my $DNA = anne_genomics::get_dna_from_fasta("$sessiondir/$queryname.fna") ;
	my %sixframes = bagel4_functions::six_frames($DNA) ;
	my @lines ;
	foreach my $frame (sort keys %sixframes) {
		push @lines, ">frame_$frame";
		push @lines, $sixframes{$frame} ;
	}	
	anne_files::write_lines("$sessiondir/$queryname.six_frames.fna", @lines) ;
	my $command = "formatdb -p T -o T -i $sessiondir/$queryname.six_frames.fna -l $sessiondir/formatdb.log"; system($command) ;
	
# 2. blast the 6 frames
	# Bactericin db I
	my ($db, $output) = ("$program_dir/$conf{bacteriocinI_proteins}", "$sessiondir/$queryname.six_frames.blast.I") ;
	$command = "blastp -outfmt 6 -query $db -db $sessiondir/$queryname.six_frames.fna -max_target_seqs 1 -num_threads $conf{cpu} -evalue $conf{blast_evalue_bacteriocinI} -out $output" ;
	system($command) ;
	&anne_files::add_header($output, "Query\tSubject\tpercent\talign_len\tmistmatches\tgap\tq_start\tq_end\ts_start\ts_end\tevalue\tbit_score");
	my %blast = anne_files::read_table_to_hash($output) ;
	
	# Bacteriocin db II
	($db, $output) = ("$program_dir/$conf{bacteriocinII_proteins}", "$sessiondir/$queryname.six_frames.blast.II") ;
	$command = "blastp -outfmt 6 -query $db -db $sessiondir/$queryname.six_frames.fna -max_target_seqs 1 -num_threads $conf{cpu} -evalue $conf{blast_evalue_bacteriocinII} -out $output" ;
	system($command) ;
	&anne_files::add_header($output, "Query\tSubject\tpercent\talign_len\tmistmatches\tgap\tq_start\tq_end\ts_start\ts_end\tevalue\tbit_score");
	%blast = (%blast, anne_files::read_table_to_hash($output)) ;
	

# 3. Get the ORFs
	my @results = "Subject\tframe\tstart\tend\tstrand\tprotein\tgene" ;
	foreach my $key (sort keys  %blast) {
		my $start = $blast{$key}{s_start} ;
		my $end = $blast{$key}{s_end} ;
		my $frame_key  ;
		if ($blast{$key}{Subject} =~ m/frame_(\d+)/) { $frame_key = $1; }  # get the frame number

		#my $prot = substr $sixframes{$frame_key},$start,$end-$start ;
		my $i;
		# extend the protein upstream until x is found
		for ($i=$start; $i>0; $i--) { last if (substr($sixframes{$frame_key}, $i, 1) eq 'x') ; }	
		$start = $i+1 ;
		# extend the protein downstream until x is found
		for ($i=$end; $i<length($sixframes{$frame_key}); $i++) { last if (substr($sixframes{$frame_key}, $i, 1) eq 'x') ; }	
		$end = $i;
		my $prot = substr $sixframes{$frame_key},$start,$end-$start ;
		print "$key\t$blast{$key}{Subject}\t$start\t$end\t$prot\n";
		
		# convert the protein position to the DNA position
		my $dna_start = $start*3 ;
		my $dna_end = $end*3 ;
		my $strand = '+' ;
		if ($blast{$key}{Subject} =~ m/frame_1/) { $dna_start +=0 ; $dna_end +=0 ; }  # correction for frame start
		if ($blast{$key}{Subject} =~ m/frame_2/) { $dna_start +=1 ; $dna_end +=1 ; }
		if ($blast{$key}{Subject} =~ m/frame_3/) { $dna_start +=2 ; $dna_end +=2 ; }
		if ($blast{$key}{Subject} =~ m/frame_(4|5|6)/) {
			$dna_start = length($DNA) - $dna_start ; 
			$dna_end   = length($DNA) - $dna_end  ; 
			$strand = '-' ;
			my $tmp = $dna_start ; $dna_start = $dna_end; $dna_end=$tmp ;  # swap start / end
			if ($1 eq '4') { $dna_start -= 0  ; $dna_end -= 0 ; }
			if ($1 eq '5') { $dna_start -= 1  ; $dna_end -= 1 ; }
			if ($1 eq '6') { $dna_start -= 2  ; $dna_end -= 2 ; }
			
		}
		
		my $gene = substr($DNA, $dna_start, $dna_end-$dna_start) ;	
		print "$dna_start\t$dna_end\t$strand\n";
		if ($strand eq '-') { $gene = anne_genomics::inverse_complement($gene); }
		my $translation = anne_genomics::translate($gene);	
		#print "Translation:\n$translation\n$prot\n" ;
		#print "$gene\n\n";
		push @results, join "\t", $key,($blast{$key}{Subject},$dna_start, $dna_end, $strand, $prot, $gene) ; 
	}
	anne_files::write_lines("$sessiondir/$queryname.blastORF", @results) ;
	


sub parseparam {
    my $var ;
    my @arg = @ARGV ;
    while(@arg) {
        $var = shift(@arg) ;
        die $usage if ($var eq '-h' or $var eq '--help') ;
		$sessiondir		= shift(@arg) if($var eq '-s') ;
		$queryname		= shift(@arg) if($var eq '-queryname') ;
    }
    die $usage if (!$queryname) ;
}


