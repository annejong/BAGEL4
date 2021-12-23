#!/usr/bin/env perl


# BAGEL4 routine to combine AOI's found by blast and primary HMMs

# 2 files will be exported: 
# - AOI.table
# - AOI.fna


use strict ;
use warnings ;
use lib "/usr/bagel4/lib" ;
use bagel4_functions ;
use lib "/usr/molgentools/lib";
use anne_files ;
use anne_genomics ;
use File::Basename;


# ---------------------------------------------------------- parameters -------------------------------------------------------------------------

my $sessiondir = '/usr/bagel4/test';
my $program_dir = dirname($0) ;
my $fna_file ;
my $queryname = 'query';
my $usage = "option:
		-s		sessiondir [default=$sessiondir]
		--dna	fasta dna sequence 
				
		this routine will combine overlapping AOIs from both AOI.blast.table and /AOI.hmm.table
		# 2 files will be exported: 
		# - AOI.table
		# - AOI.fna
				
		e.g.  /usr/bagel4/bagel4_combine_regions.pl -s /usr/bagel4/test_plantarum -dna /var/genomes/g2d_mirror/Streptococcus_pneumoniae_R6/ASM704v1_genomic.fna	
		/usr/bagel4/bagel4_combine_regions.pl -s /usr/bagel4/test -dna /usr/bagel4/test/NC_0085331.fna -queryname NC_0085331
" ;


&parseparam() ;

# ---------------------------------------------------------- main -------------------------------------------------------------------------
print "---------------------------------------------- Combine overlapping AOIs ------------------------------------------------------------\n";

my @lines = anne_files::read_lines("$sessiondir/$queryname.AOI.blast.table") ;
push @lines, anne_files::read_lines("$sessiondir/$queryname.AOI.hmm.table") ;

my %table ;
my $key = 0;
print "\t Read AOIs from hmm and blast identification\n";
foreach my $line (@lines) {
	my @items = split /\t/, $line ;
	
	if (scalar @items > 4) {
		if ($items[1] > $items[2]) { my $tmp=$items[1]; $items[1]=$items[2]; $items[2]=$tmp; }  # swap if start < end
		$key++ ;
		$table{$key}{start} 	= $items[1] ;
		$table{$key}{end} 		= $items[2] ;
		$table{$key}{mid} 		= $table{$key}{start} + ($table{$key}{end}-$table{$key}{start})/2 ;
		$table{$key}{strand}	= $items[3] ;
		$table{$key}{name} 		= $items[4] ;
		$table{$key}{type} 		= $items[5] ;
		$table{$key}{bit_score} = 0 ;
		$table{$key}{bit_score} = $items[6] if (defined($items[6])) ;
		$table{$key}{line} 		= $line ;
		my $size = $table{$key}{end}  - $table{$key}{start} ;
		print "$table{$key}{name}\t$table{$key}{start}\t$table{$key}{end}\t$table{$key}{strand}\tSize=$size\n";
		
	}	
}

my @delkeys ;
my @sorted_keys ;
foreach my $key (sort {$table{$a}{start} <=> $table{$b}{start}} keys %table) {
	push @sorted_keys, $key ;
}

for (my $i = 0; $i<(scalar @sorted_keys)-1; $i++) {
	print "=========>".$sorted_keys[$i]."\n" ;
	if ($table{$sorted_keys[$i+1]}{mid} <= $table{$sorted_keys[$i]}{end}) {
		push @delkeys, $sorted_keys[$i] ;
		$table{$sorted_keys[$i+1]}{start} = $table{$sorted_keys[$i]}{start} ;
		if ($table{$sorted_keys[$i]}{bit_score} > $table{$sorted_keys[$i+1]}{bit_score} ) {
			$table{$sorted_keys[$i+1]}{name} = $table{$sorted_keys[$i]}{name} ;
			$table{$sorted_keys[$i+1]}{bit_score} = $table{$sorted_keys[$i]}{bit_score} ;
		} 
	}
}


# OLD
# for (my $iteration = 1; $iteration<5; $iteration++) {  
# 	for (my $i = 1; $i<$n; $i++) { # compare each $i (1..n-1) against all $j (n+1..n)
# 		for (my $j = $i+1; $j<=$n; $j++) {  
# 			if ($table{$j}{start}>=$table{$i}{start} and $table{$j}{start}<=$table{$i}{end}) {  # adjust end position of $i and remove $j
# 				$table{$i}{end} = $table{$j}{end}  ;
# 				$table{$i}{name} = $table{$j}{name} if ($table{$i}{bit_score} < $table{$j}{bit_score}) ;
# 				push @delkeys, $j ;
# 				print "\t\toverlap between $i and $j\n";
# 			}	
# 			if ($table{$j}{end}>=$table{$i}{start} and $table{$j}{end}<=$table{$i}{end}) { # adjust start position of $i
# 				$table{$i}{start} = $table{$j}{start} ;
# 				$table{$i}{name} = $table{$j}{name} if ($table{$i}{bit_score} < $table{$j}{bit_score})  ; 
# 				push @delkeys, $j ;
# 				print "\t\toverlap between $i and $j\n";
# 			}	
# 		}		
# 	}
# }


@delkeys = anne_files::unique_array(@delkeys);
print "\nNumber of AOIs merged:".(scalar @delkeys)."\n";

my $DNA = anne_genomics::get_dna_from_fasta($fna_file) ;
my $DNA_len = length($DNA);
my @AOI_table ;
my @fna ;
my $AOI_number = 0 ;
foreach my $key (sort keys %table) {
	 if ( !grep( /^$key$/, @delkeys ) ) { 
		$AOI_number++ ;
		push @fna, ">$queryname.AOI_0".$AOI_number ;
		push @fna, substr $DNA, $table{$key}{start} , $table{$key}{end} -$table{$key}{start} ;
		my @row = "$queryname.AOI_0".$AOI_number ;
		push @row, $DNA_len ;
		push @row, $table{$key}{start} ;
		push @row, $table{$key}{end} ;
		push @row, $table{$key}{strand};
		$table{$key}{name} = anne_files::unique_pipe_list($table{$key}{name}) ;
		push @row, $table{$key}{name} ;	
		push @row, $table{$key}{end} - $table{$key}{start} ;	
		print join("\t", @row) ;
		print "\n";
		push @AOI_table, join("\t", @row) ;
	}	
}	 
anne_files::write_lines("$sessiondir/$queryname.AOI.table", @AOI_table);
anne_files::write_lines("$sessiondir/$queryname.AOI.fna", @fna);
# anne_files::table2html("$sessiondir/$queryname.AOI.table", "$sessiondir/$queryname.AOI.table.html");

print "Results written to table $queryname.AOI.table\n";
print "Results written to fasta $queryname.AOI.fna\n";

# ---------------------------------------------------------- functions -------------------------------------------------------------------------


sub parseparam {
    my $var ;
    my @arg = @ARGV ;
    while(@arg) {
        $var = shift(@arg) ;
        die $usage if ($var eq '-h' or $var eq '--help') ;
		$sessiondir		= shift(@arg) if($var eq '-s') ;
		$fna_file		= shift(@arg) if($var eq '-dna') ;
		$queryname		= shift(@arg) if($var eq '-queryname') ;
    }
	die $usage if (!$fna_file) ;
}


