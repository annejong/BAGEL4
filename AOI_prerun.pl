#!/usr/bin/env perl

use strict ;
use warnings ;
use File::Basename;
use lib "/usr/molgentools/lib";
use anne_files ;
use anne_genomics ;
use anne_misc;
use lib "/data/bagel4/lib" ;
use bagel4_functions ;


my $usage = "
This routine is made to speed-up WGS or metagenome data containing large number reads above the cutoff size but relatively short (typically 2000 - 20.000 bp)
By quick initial screening for proximity HMM and a multi record blast, enables filtering of reads without a potential AOI

example for tesing: 
cd /tmpdrive/bagel4wrapper/129.125.142.95d2jlafq9ugmgpdb54ddpllbna7784
sudo /data/bagel4/AOI_prerun.pl -s /tmpdrive/bagel4wrapper/129.125.142.95d2jlafq9ugmgpdb54ddpllbna7784 -dna GCF_000177595.1_ASM17759v1.fna
sudo /data/bagel4/AOI_prerun.pl -s /tmp/BAGEL4WRAPPER/129.125.142.95ucgs78u8e2r8eo13f0ehivvl80871 -dna GCF_001123405.1_6259_5_20.fna
";

# ---------------------------------------------------------- parameters -------------------------------------------------------------------------

my $program_dir = dirname($0) ;
my $sessiondir = '.';
my $fna_file ;

parseparam();

# --------------------------------------------------------------------- main -------------------------------------------------------------------------------	

my %conf    = bagel4_functions::read_conf("$program_dir/bagel4.conf") ;
my %contigs = anne_genomics::fasta2hash("$sessiondir/$fna_file");
my $six_frames_multifasta_filename = "$sessiondir/01.prescan_proteindb.faa" ;


# 1. Translate all records in the fasta file
	translate_multi_fasta() ;


# 2. BLAST all 6 protein frames
	my $command = "blastp -outfmt 6 -db $program_dir/$conf{bacteriocinI_proteins} -query $six_frames_multifasta_filename -max_target_seqs 1 -num_threads $conf{cpu} -evalue $conf{blast_evalue_bacteriocinI}   -out $sessiondir/01.prescan.blast1" ;
	print "===>$command\n\n"; system($command) ;
	$command = "blastp -outfmt 6 -db $program_dir/$conf{bacteriocinII_proteins}   -query $six_frames_multifasta_filename -max_target_seqs 1 -num_threads $conf{cpu} -evalue $conf{blast_evalue_bacteriocinII}  -out $sessiondir/01.prescan.blast2" ;
	print "===>$command\n\n"; system($command) ;
	$command = "blastp -outfmt 6 -db $program_dir/$conf{bacteriocinIII_proteins}  -query $six_frames_multifasta_filename -max_target_seqs 1 -num_threads $conf{cpu} -evalue $conf{blast_evalue_bacteriocinIII} -out $sessiondir/01.prescan.blast3" ;
	print "===>$command\n\n"; system($command) ;
	$command = "cat $sessiondir/01.prescan.blast1 $sessiondir/01.prescan.blast2 $sessiondir/01.prescan.blast3 > $sessiondir/01.prescan.blast123";	
	system($command);
	
# 3. proximity HMM search in all 6 frames. Make one hmms file of all proximity singles. locate in subfolder /db_hmm
		# cat CyaG.hmm       >all_proximity.hmm
		# cat labKC.hmm     >>all_proximity.hmm
		# cat LinL.hmm      >>all_proximity.hmm
		# cat PF00881.hmm   >>all_proximity.hmm
		# cat PF01944.hmm   >>all_proximity.hmm
		# cat PF02794.hmm   >>all_proximity.hmm
		# cat PF04055.hmm   >>all_proximity.hmm
		# cat PF05147.hmm   >>all_proximity.hmm
		# cat PF13471.hmm   >>all_proximity.hmm
		# cat TIGR04195.hmm >>all_proximity.hmm

	my $multi_HMMs_file = "$program_dir/db_hmm/all_proximity.hmm";
	my $domtblout = "01.prescan_proximity.domtblout";
	$command = "$conf{hmmsearch}/hmmsearch --noali --cpu $conf{cpu} --domtblout $sessiondir/$domtblout --domT $conf{domT} $multi_HMMs_file $six_frames_multifasta_filename" ;
	print "===>$command\n\n"; system($command) ;


# 4. remove fasta records without any hit
	my @blasthits = anne_files::read_lines("$sessiondir/01.prescan.blast123") ;	
	my @HMMhits = anne_files::read_lines("$sessiondir/$domtblout") ;	
	my @to_keep ;
	foreach my $line (@blasthits) {	if ($line =~ m/^(.+?)\t/) { push @to_keep, $1; } }
	foreach my $line (@HMMhits) { if ($line !~ m/^#/ and $line =~ m/^(.+?)\s+/) { push @to_keep, $1; } }
	@to_keep = anne_files::unique_array(@to_keep) ;
	print "Record to keep =".(scalar @to_keep)."\n"; ;
	foreach my $key (@to_keep) { print "to keep:  $key\n";	}
	my @lines ;
	foreach my $contig (keys %contigs) {
		if (anne_misc::str_in_array($contig, @to_keep)) {
			push @lines, ">$contig\n$contigs{$contig}\n";
		}
	}
	anne_files::write_lines($fna_file, @lines) ;
	
# --------------------------------------------------------------------- functions ----------------------------------------------------------------------------	
	
sub translate_multi_fasta {
	
	my @lines ;
	foreach my $contig (keys %contigs) {
		my %sixframes = bagel4_functions::six_frames($contigs{$contig}) ;
		foreach my $frame (keys %sixframes) {
			#push @lines, '>'.$contig.'_frame'.$frame."\n$sixframes{$frame}\n";
			push @lines, ">$contig\n$sixframes{$frame}\n";
		}
	}
	anne_files::write_lines($six_frames_multifasta_filename, @lines) ;
	
}	

sub parseparam {
    my $var ;
    my @arg = @ARGV ;
    while(@arg) {
        $var = shift(@arg) ;
        die $usage if ($var eq '-h' or $var eq '--help') ;
		$sessiondir		= shift(@arg) if($var eq '-s') ;
		$fna_file		= shift(@arg) if($var eq '-dna') ;
    }
    die $usage if (!$fna_file) ;
}

	