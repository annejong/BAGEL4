#!/usr/bin/env perl

# convert GeneTable to Genbank


use strict ;
use warnings ;
use lib "/data/bagel4/lib" ;
use bagel4_functions ;
use lib "/data/molgentools/lib";
use anne_files ;
use File::Basename;

# ---------------------------------------------------------- parameters -------------------------------------------------------------------------

my $sessiondir = '.';
my $program_dir = dirname($0) ;
my $queryname  ;
my $outputfilename = 'my_genetable.gbk';


my $usage = "option:
		-i		GeneTable
		-o		output filename
				
		convert GeneTable to Genbank
				
		e.g.  /data/bagel4/bagel4_GeneTable_2_genbank.pl -i NC_0085331.AOI_01
" ;


&parseparam() ;
my $indent = "                     ";
# read the GeneTable
#region_name<--->region_size<--->orf<--->gene_name<----->gene_start<---->gene_end<------>gene_strand<--->real_start<---->real_end<------>real_strand<--->gene_color<---->function<------>motifs<>annotation<---->protein>dna
my %table = anne_files::Table2hash_v2("$sessiondir/$queryname.GeneTable", 'true') ;


anne_files::write_lines("$sessiondir/$queryname.gbk",GeneTable_2_Genbank() ) ;

#print "\t =====> alignment file for $query written to $outputfilename\n";



# ---------------------------------------------------------- main -------------------------------------------------------------------------------

sub format_dna {  # make the genbank DNA format
	my $dna = shift ;
	my @results ;
	for (my $i = 1; $i <= length($dna); $i +=50) {
		my $row = substr $dna, $i, 50 ;
		$row =~ s/..........\K(?=.)/ /sg;
		push @results, sprintf ("%9s", $i).' '.$row ;
	}
	my $result = join("\n", @results) ;
	return $result;
}	

sub format_feature {
	# Comply to the nasty Genbank layout
	my $prot = shift ;  # including feature name e.g. /translation=
	my @results ;
	for (my $i = 0; $i <= length($prot); $i +=58) {
		my $row = substr $prot, $i, 58 ;
		push @results, $indent.$row ;
	}
	return join("\n", @results) ;
}

sub GeneTable_2_Genbank {
	
	my @genbank = "LOCUS       $table{1}{region_name}            $table{1}{region_size} bp    DNA     linear";
	push @genbank, "FEATURES             Location/Qualifiers
     source          1..$table{1}{region_size}
                     /organism=\"unkown\"
                     /mol_type=\"genomic DNA\"
                     /strain=\"unkown\"" ;
	foreach my $key (sort {$table{$a}{gene_start} <=> $table{$a}{gene_start}} keys %table) {
		push @genbank, "     gene            $table{$key}{gene_start}..$table{$key}{gene_end}"  if ($table{$key}{gene_strand} eq '+') ;
		push @genbank, "     gene            complement($table{$key}{gene_start}..$table{$key}{gene_end})"  if ($table{$key}{gene_strand} eq '-') ;
		push @genbank, "$indent/gene=\"$table{$key}{gene_name}\"" ;
		push @genbank, "$indent/locus_tag=\"$table{$key}{orf}\"" ;
		push @genbank, "     CDS             $table{$key}{gene_start}..$table{$key}{gene_end}"  if ($table{$key}{gene_strand} eq '+') ;
		push @genbank, "     CDS             complement($table{$key}{gene_start}..$table{$key}{gene_end})"  if ($table{$key}{gene_strand} eq '-') ;
		push @genbank, "$indent/gene=\"$table{$key}{gene_name}\"" ;
		push @genbank, "$indent/locus_tag=\"$table{$key}{orf}\"" ;
		push @genbank, format_feature('/product="'.$table{$key}{function}.'"') ;
		push @genbank, format_feature('/note="'.$table{$key}{annotation}.'"') ;
		push @genbank, format_feature('/translation="'.$table{$key}{protein}.'"') ;
	}
	push @genbank, "ORIGIN" ;
	my $dna = anne_genomics::get_dna_from_fasta("$sessiondir/$queryname.fna") ;
	push @genbank, format_dna($dna);
	push @genbank,'//';
	return @genbank ;

}

# ------------------------------------------------------------ functions ---------------------------------------------------------------------------


sub parseparam {
    my $var ;
    my @arg = @ARGV ;
    while(@arg) {
        $var = shift(@arg) ;
        die $usage if ($var eq '-h' or $var eq '--help') ;
		$sessiondir	= shift(@arg) if($var eq '-s') ;
		$queryname	= shift(@arg) if($var eq '-i') ;
		$outputfilename	= shift(@arg) if($var eq '-o') ;
    }
    die $usage if (!$queryname) ;
}

