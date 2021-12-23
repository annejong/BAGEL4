#!/usr/bin/env perl

# After identification of the AOIs on the basis of HMM and BLAST, this routine will annotate the AOI


use strict ;
use warnings ;
use lib "/data/bagel4/lib" ;
use bagel4_functions ;
use lib "/data/molgentools/lib";
use anne_files ;
use anne_genomics ;
use File::Basename;


# ---------------------------------------------------------- parameters -------------------------------------------------------------------------

my $sessiondir = '/data/bagel4/test';
my $program_dir = dirname($0) ;
my $cpu = 16 ;
my $fna_file = 'AOI.fna' ;
my $queryname = 'query' ;


my %conf = bagel4_functions::read_conf("$program_dir/bagel4.conf") ;

my $usage = "option:
		-s		sessiondir [default=$sessiondir]
		-dna	dna all AOIs in FASTA format
		-cpu	number of cpu's used [default=$cpu]
		-query  name of the query
				
		this routine will perform the annotation of the identified AOI
				
		e.g.  time /usr/bagel4/bagel4_annotation_of_AOI.pl -s /usr/bagel4/test -model false -dna /usr/bagel4/test/NC_0085331.AOI_5.fna -query NC_0085331
		e.g.  time /usr/bagel4/bagel4_annotation_of_AOI.pl -s /usr/bagel4/test_plantarum -model true -dna AOI.fna
" ;


&parseparam() ;


my $href_bagel4 = 'http://bagel4.molgenrug.nl/bagel4results';
my $tmp_folder = '/tmp/BAGEL4WRAPPER' ;
my $sessionID = $sessiondir ;
$sessionID =~ s/$tmp_folder\/// ;


# ---------------------------------------------------------- main -------------------------------------------------------------------------------



# 1. load descriptions of pfam and bagel4_hmm models
	my %hmm_descriptions = anne_files::read_table_to_hash("$program_dir/$conf{bacteriocin_hmm_db_descriptions}") ; # load bacteriocin_hmm_db_descriptions
							
# 2. Annotation of the AOI's	
	my %GeneTable ;  # container for all the final gene table 
	my %fasta = anne_genomics::fasta2hash($fna_file) ;
	my %AOI_sizes;
	# the color rules are needed to color the context genes
	my %color_rules = anne_files::read_table_to_hash("$program_dir/tables/bagel4_context_annotation_and_color_rules.txt"); # header: key	rule	min_len	max_len	name	color
	my %swissport_descriptions = anne_files::read_table_to_hash("$program_dir/swissprot/bagel4_swissport_descriptions_pfams.table"); # header: UniRefID	Function	Tax	PFAMS
	my %AOI_table = read_AOI_table("$sessiondir/$queryname.AOI.table") ;

	foreach my $key (sort keys %fasta) {
		# ORF calling
			print "Annotation of $key\n";
			bagel4_functions::write_bagel_log( "ORF calling on AOI $key", $sessiondir,1,1,1 );
			anne_files::write_string("$sessiondir/$key.fna" , ">$key\n$fasta{$key}\n") ;
			# the fasta DNA file
			call_genes($key, "$sessiondir/$key.fna" ) ;	# use glimmer to call genes and glimmer model
			$AOI_sizes{$key} = length($fasta{$key}) ;  # length of the AOI
			glimmer_2_GeneTable($key) ; 	# add the genes predicted by Glimmer to the GeneTable
			call_small_orfs($key) ;			# NOTE: NOT for PKS NRPS Terpene and BacteriocinIII
			translate_genes($key) ;			# translate the genes to proteins and Add to GeneTable
			add_original_blast($key) ;		# if the ORF callers fail to find the original BLAST hit on which the AOI was found, this routine recover this
		# Annotation
			bacteriocin_hmmsearch($key) ;		# Functional annotation on the basis of the bacteriocin_hmm's  <============================================================= check later 
			bacteriocin_hmm_2_GeneTable($key, $conf{bacteriocinhmmColor}, "$sessiondir/$key.bacteriocin_hmmsearch"); 		# data to the GeneTable   <=========================== check later 
			blast_bacteriocin_db($key) ;			# blast the all orfs against the bacteriocin database and add to GeneTable
			#remove_unannotated_sORFS($key) ;
			blast_SwissProt($key, "$sessiondir/$key.faa", "$sessiondir/$key.blast_RefSeq90.table") ;		# blast the proteins against the SwissProt UniRef90 database
			SwissProt_2_GeneTable($key,  "$sessiondir/$key.blast_RefSeq90.table" ) ;
			add_real_coords_2_GeneTable($key) ;
			save_GeneTable($key, "$sessiondir/$key.GeneTable") ;				# save the genetable of the key
			# save_GeneTable_asHTML($key) ;		# save the genetable of the key as html
			# save_Alignments_asHTML($key, "$sessiondir/$key.Alignments.html" ) ;		# save all alignment files as html
			my $command = "$program_dir/bagel4_GeneTable_2_genbank.pl -s $sessiondir -i $key" ; system($command) ; 
		# adding RNA_seq data
			
				
			
	}
	GeneTable_2_ptt() ; # write a ptt file for all $key's
	save_all_GeneTables("$sessiondir/$queryname.GeneTables");
	# save_all_GeneTables_asHTML("$sessiondir/$queryname.GeneTables.html");
	# save_all_GeneTables_asHTML("$sessiondir/00.GeneTables.html");  # For multiple FNA input this should be updated

	
# ------------------------------------------------------------ functions ---------------------------------------------------------------------------

sub add_original_blast {
	# Blast db against 6 frames of the AOI and add the protein if not present in de GeneTable
	# BlastTable = "Subject\tframe\tstart\tend\tstrand\tprotein\tgene" ;
	my $key = shift ;
	my $command = "$program_dir/bagel4_ORF_from_frameblast.pl -s $sessiondir -queryname $key" ;
	system($command) ;
	my %blast = anne_files::read_table_to_hash("$sessiondir/$key.blastORF") ;
	my @fasta =  anne_files::read_lines("$sessiondir/$key.faa") ;
	my $count = 0;
	foreach my $blast_key (sort keys %blast) {
		print "=====> add_original_blast: $key\n";
		# gene end not in GeneTable, add this to the genetabl
		my $match = 0 ;
		foreach my $orf (sort keys %{$GeneTable{$key}}) {
			if ( $blast{$blast_key}{strand} eq '+' and $GeneTable{$key}{$orf}{gene_end} > $blast{$blast_key}{end}-5     and $GeneTable{$key}{$orf}{gene_end} < $blast{$blast_key}{end}+5 )     { $match = 1; last; }
			if ( $blast{$blast_key}{strand} eq '-' and $GeneTable{$key}{$orf}{gene_start} > $blast{$blast_key}{start}-5 and $GeneTable{$key}{$orf}{gene_start} < $blast{$blast_key}{start}+5 ) { $match = 1; last; }
		}
		if (!$match) {
			$count++ ;
			$GeneTable{$key}{$blast_key}{region_name} = $key ;
			$GeneTable{$key}{$blast_key}{region_size} = $AOI_sizes{$key} ;
			$GeneTable{$key}{$blast_key}{orf}	 		= 'orfblast_'.$count ;
			$GeneTable{$key}{$blast_key}{gene_name} 	= $blast_key ;
			$GeneTable{$key}{$blast_key}{gene_start} 	= $blast{$blast_key}{start} ;
			$GeneTable{$key}{$blast_key}{gene_end} 		= $blast{$blast_key}{end};
			$GeneTable{$key}{$blast_key}{gene_strand}	= $blast{$blast_key}{strand} ;
			$GeneTable{$key}{$blast_key}{gene_color} 	= $conf{bacteriocinIColor} ;
			$GeneTable{$key}{$blast_key}{function} 		= $blast_key ;
			$GeneTable{$key}{$blast_key}{motifs} 		= '' ;
			$GeneTable{$key}{$blast_key}{alignment_url} = '';
			$GeneTable{$key}{$blast_key}{alignment_filename} = '';
			$GeneTable{$key}{$blast_key}{annotation} 	= '' ;
			$GeneTable{$key}{$blast_key}{dna}         	= $blast{$blast_key}{gene};
			$GeneTable{$key}{$blast_key}{protein} 		= $blast{$blast_key}{protein};
			$GeneTable{$key}{$blast_key}{comment} 		.= 'database hit';
			# add protein to the fast file
			push @fasta, ">orfblast_$count\n$blast{$blast_key}{protein}" ;
		}
	}
	anne_files::write_lines("$sessiondir/$key.faa", @fasta); # Write the all proteins including the added proteins
}


sub remove_unannotated_sORFS {
	print "---------------------------------------- remove unannotated sORFs ---------------------------\n";
	my $key = shift ;
	## my $count = 0 ;
	## my @colors = ($conf{bacteriocinIColor},$conf{bacteriocinIIColor},$conf{bacteriocinIIIColor},$conf{bacteriocinhmmColor}) ;
	## foreach my $orf (sort keys %{$GeneTable{$key}}) { $count++ if (anne_files::in_array($GeneTable{$key}{$orf}{gene_color},@colors)) ; }
	## # if at least one green colored gene is found we remove the sORFs without function
	foreach my $orf (sort keys %{$GeneTable{$key}}) { 
		if ($orf =~ m/^sORF/ and $GeneTable{$key}{$orf}{function} eq '') { delete $GeneTable{$key}{$orf}; }
	}
	# remove the ORFs from the protein file
	my @faa ;
	foreach my $orf (sort keys %{$GeneTable{$key}}) { push @faa, ">$orf\n$GeneTable{$key}{$orf}{protein}"; }
	anne_files::write_lines("$sessiondir/$key.faa",@faa );
}

sub save_Alignments_asHTML {
	# save as HTML table
	my ($key,$outfilename) = @_ ;
	my @html  ;
	foreach my $orf (sort keys %{$GeneTable{$key}}) {
		push @html, anne_files::read_lines($GeneTable{$key}{$orf}{alignment_filename}) if ($GeneTable{$key}{$orf}{alignment_filename} ne '') ;	
	}
	anne_files::write_lines($outfilename, @html);
}



sub read_AOI_table {
	my @lines = anne_files::read_lines(shift) ;
	my %results ;
	foreach my $line (@lines) {
		my @items = split /\t/, $line ;
		if (scalar @items > 5) {
			$results{$items[0]}{start} = $items[2] ;
			$results{$items[0]}{end} = $items[3] ;
			$results{$items[0]}{strand} = $items[4] ;
			$results{$items[0]}{size} = $items[6] ;
		}
	}	
	return %results ;
}

sub add_real_coords_2_GeneTable {
	my $key = shift ;
	foreach my $orf (sort keys %{$GeneTable{$key}}) {
		if ($AOI_table{$key}{strand} eq '+') {
			$GeneTable{$key}{$orf}{real_start} 	= $GeneTable{$key}{$orf}{gene_start} + $AOI_table{$key}{start}  ;
			$GeneTable{$key}{$orf}{real_end} 	= $GeneTable{$key}{$orf}{gene_end}   + $AOI_table{$key}{start}  ;
			$GeneTable{$key}{$orf}{real_strand} = $GeneTable{$key}{$orf}{gene_strand} ;
		} else {
			$GeneTable{$key}{$orf}{real_start} 	= $AOI_table{$key}{end} - $GeneTable{$key}{$orf}{gene_end} ;
			$GeneTable{$key}{$orf}{real_end} 	= $AOI_table{$key}{end} - $GeneTable{$key}{$orf}{gene_start} ;
			$GeneTable{$key}{$orf}{real_strand} = '-' if ( $GeneTable{$key}{$orf}{gene_strand} eq '+') ;  # swap the real gene strand
			$GeneTable{$key}{$orf}{real_strand} = '+' if ( $GeneTable{$key}{$orf}{gene_strand} eq '-') ;
		}
	}	
}	
	


sub call_small_orfs {
	my $key = shift ;
	my $command = "$program_dir/bagel4_small_orf_calling.pl -s $sessiondir -queryname $key  -min $conf{smallorf_min_prot_length} -max $conf{smallorf_max_prot_length}";
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
	print "$command\n";
	print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n";
	system($command) ;
	my %small_orfs = anne_files::read_table_to_hash("$sessiondir/$key.sOFS.table") ;
	#foreach my $orf (sort {$GeneTable{$a}{start} <=> $GeneTable{$b}{start}} keys %small_orfs) {
	foreach my $orf (keys %small_orfs) {
		$GeneTable{$key}{$orf}{region_name} 	= $key ;
		$GeneTable{$key}{$orf}{region_size} 	= $AOI_sizes{$key} ;
		$GeneTable{$key}{$orf}{orf}	 		= $orf ;
		$GeneTable{$key}{$orf}{gene_name} 	= $orf ;
		$GeneTable{$key}{$orf}{gene_start} 	= $small_orfs{$orf}{start} ;
		$GeneTable{$key}{$orf}{gene_end} 	= $small_orfs{$orf}{end} ;
		$GeneTable{$key}{$orf}{gene_strand}	= $small_orfs{$orf}{strand} ;
		$GeneTable{$key}{$orf}{gene_color} 	= $conf{defaultColor} ;
		$GeneTable{$key}{$orf}{function} 	= '' ;
		$GeneTable{$key}{$orf}{motifs} 		= "RBS=$small_orfs{$orf}{rbs}<br>"  ;
		$GeneTable{$key}{$orf}{alignment_url} = '';
		$GeneTable{$key}{$orf}{alignment_filename} = '';
		$GeneTable{$key}{$orf}{annotation} 	= '' ;
		$GeneTable{$key}{$orf}{dna}         = $small_orfs{$orf}{gene} ;
		$GeneTable{$key}{$orf}{protein} 	= anne_genomics::translate($small_orfs{$orf}{gene});
		$GeneTable{$key}{$orf}{comment} 	.= 'small ORF';
	}
	# add the small genes to all genes
	my @lines = (anne_files::read_lines("$sessiondir/$key.genes"),anne_files::read_lines("$sessiondir/$key.sOFS.fna")) ;
	anne_files::write_lines("$sessiondir/$key.genes", @lines);
}

sub save_GeneTable {
	# save as tab delimited table
	my ($key,$outfilename) = @_ ;
	my @headers = ('region_name','region_size','orf','gene_name','gene_start','gene_end','gene_strand','real_start','real_end','real_strand','gene_color','function','motifs','annotation','protein','dna');
	my @lines = join "\t", @headers ;
	foreach my $orf (sort keys %{$GeneTable{$key}}) {
		my @row ;
		foreach my $header (@headers) { push @row,$GeneTable{$key}{$orf}{$header} ; }
		push @lines, join "\t", @row ;
	}
	anne_files::write_lines($outfilename, @lines);
}


sub save_GeneTable_asHTML {
	# save as HTML table
	my $key = shift ;
	# the full table
	my @headers = ('ORF','gene_name','start','end','strand','real_start','real_end','real_strand','function','motif','annotation','protein');
	my @html = "<a id=GeneTable onclick=toggleTable(\'$key\'); href=javascript:void(0);>Full GeneTable $key</a><br>" ;
	push @html, "<TABLE id=$key class=display><tr><th bgcolor=black;color=white>".(join "</th><th>", @headers).'</th></tr>' ;
	# the condensed table
	my @cond_headers = ('Gene','Function','Motif');
	my @cond_html; 
	#= "<a id=GeneTableSmall onclick=toggleTable(\'cond_$key\'); href=javascript:void(0);>GeneTable $key</a><br>" ;
	push @cond_html, "<TABLE id=cond_$key class=display><tr><th bgcolor=black;color=white>".(join "</th><th>", @cond_headers).'</th></tr>' ;
	foreach my $orf (sort {$GeneTable{$key}{$a}{gene_start} <=> $GeneTable{$key}{$b}{gene_start}} keys %{$GeneTable{$key}}) {
		# the full table
		$GeneTable{$key}{$orf}{annotation} .=  $GeneTable{$key}{$orf}{alignment_url} ;	
		my  @row = ($GeneTable{$key}{$orf}{orf},$GeneTable{$key}{$orf}{gene_name}) ;
		push @row, ($GeneTable{$key}{$orf}{gene_start},$GeneTable{$key}{$orf}{gene_end},$GeneTable{$key}{$orf}{gene_strand}) ;
		push @row, ($GeneTable{$key}{$orf}{real_start},$GeneTable{$key}{$orf}{real_end},$GeneTable{$key}{$orf}{real_strand}) ;
		push @row ,($GeneTable{$key}{$orf}{function},$GeneTable{$key}{$orf}{motifs},$GeneTable{$key}{$orf}{annotation}, $GeneTable{$key}{$orf}{protein}) ;
		push @html, "<tr align=left bgcolor=$GeneTable{$key}{$orf}{gene_color}><td>".(join "</td><td>", @row).'</td></tr>' ;
		# the condensed table
		my @cond_row = ($GeneTable{$key}{$orf}{gene_name},$GeneTable{$key}{$orf}{function},$GeneTable{$key}{$orf}{motifs}) ;
		push @cond_html,"<tr align=left bgcolor=$GeneTable{$key}{$orf}{gene_color}><td>".(join "</td><td>", @cond_row).'</td></tr>' ;
	}
	push @html, '</table>' ;
	push @cond_html, '</table>' ;
	anne_files::write_lines("$sessiondir/$key.GeneTable.html", @html);
	anne_files::write_lines("$sessiondir/$key.GeneTableCondensed.html", @cond_html);
}



sub save_all_GeneTables_asHTML {
	# make a html overview of all GeneTables_xxx.html
	my $outfilename = shift ;
	my @html = "Hide / Show Gene Table<br>" ;
	foreach my $key (sort keys %GeneTable) {
		push @html, anne_files::read_lines("$sessiondir/$key.GeneTable.html") ;
	}	
	anne_files::write_lines($outfilename, (@html));
}


sub save_all_GeneTables {
	my $outfilename = shift ;
	my @headers = ('region_name','region_size','orf','gene_name','gene_start','gene_end','gene_strand','gene_color','function','motifs','annotation');
	my @lines = join "\t", @headers ;
	foreach my $key (sort keys %GeneTable) {
		foreach my $orf (sort {$GeneTable{$key}{$a}{gene_start} <=> $GeneTable{$key}{$b}{gene_start}} keys %{$GeneTable{$key}}) {
		# foreach my $orf (sort keys %{$GeneTable{$key}}) {
			my @row ;
			foreach my $header (@headers) { push @row,$GeneTable{$key}{$orf}{$header} ; }
			push @lines, join "\t", @row ;
		}	
	}
	anne_files::write_lines($outfilename, @lines);
}



sub blast_SwissProt  {
	# blast the proteins against the SwissProt UniRef90 database
	my $key = shift ;  # the queryname
	my $proteinfile = shift ;
	my $outfile = shift ;
	my $command = "$program_dir/bagel4_blast_swissprot.pl -s $sessiondir -db $program_dir/swissprot/archea_bacteria -evalue 1E-03 $proteinfile -query $sessiondir/$key.faa -cpu $conf{cpu} -out $outfile" ;
	system($command) ;
}


sub SwissProt_2_GeneTable {
	# Use the blast_SwissProt to add the annotation
	# header: Query	UniProt	Percent	Evalue	Function	Tax	PFAM	url
	# use color rules of bagel4_context_annotation_and_color_rules.txt
	my ($key, $inputfile) = @_ ;
	my @result ;
	my %blast_results = anne_files::read_table_to_hash($inputfile) ;	# header: Query	UniProt	Percent	Evalue	Function	Tax	PFAM	url
	foreach my $orf (sort keys %blast_results) {
		if ($GeneTable{$key}{$orf}{gene_color} ne $conf{bacteriocinIColor}) {
		#$GeneTable{$key}{$orf}{gene_color} 	= get_color_from_color_rules('UniRef90_'.$blast_results{$orf}{UniProt}, $blast_results{$orf}{PFAMS}) ;  # add the UniRef90 tag
			my $ColorRuleKey = get_color_from_color_rules($blast_results{$orf}{UniProt}, $blast_results{$orf}{PFAMS}, length($GeneTable{$key}{$orf}{protein}) ) ;  # get the key from the table "bagel4_context_annotation_and_color_rules.txt"
			if ($ColorRuleKey ne '') { 
				$GeneTable{$key}{$orf}{gene_color} 	= $color_rules{$ColorRuleKey}{color} ; 
				$GeneTable{$key}{$orf}{colorrulename} 	= $color_rules{$ColorRuleKey}{name} ; 
				$GeneTable{$key}{$orf}{gene_name} 	= $color_rules{$ColorRuleKey}{name} ; 
			} else {
				$GeneTable{$key}{$orf}{gene_color} 	= $conf{UniProtColor} ;
			}
			$GeneTable{$key}{$orf}{motifs} 		.= $blast_results{$orf}{PFAMS} ;
			$GeneTable{$key}{$orf}{function} 	= $blast_results{$orf}{Function} ;
			$blast_results{$orf}{url} =~ s/UniRef90_//g ;
			$GeneTable{$key}{$orf}{annotation} 	= "Species=$blast_results{$orf}{Tax}<br>match=$blast_results{$orf}{Percent}%|Evalue=$blast_results{$orf}{Evalue}<br>UniRef:$blast_results{$orf}{url}" ;
		}
	}
}

sub get_color_from_color_rules {
	my $UniRef = shift ;  # uniprot ID or context protein ID
	my $query_pfams = shift ;
	my $len = shift ; # length of the protein
	my $ColorRuleKey = '' ;
	if ($query_pfams ne '') {
		foreach my $key (sort {$a <=> $b} keys %color_rules) { 
			my @rule_pfams = split ";", $color_rules{$key}{rule} ;
			my $count = 0 ;
			foreach my $rule_pfam (@rule_pfams) { $count++ if ($query_pfams =~ m/$rule_pfam/) ; }
			if ($count >= scalar @rule_pfams and $len>=$color_rules{$key}{min_len} and $len<=$color_rules{$key}{max_len} )   { $ColorRuleKey = $key ; }
		}
	}
	#print "===>get_color_from_color_rules==> $UniRef;$query_pfams;colorkey=$ColorRuleKey;len=$len \n";
	return $ColorRuleKey ;
}

sub read_glimmer_predict {
	my $key = shift ;
	my @lines = anne_files::read_lines("$sessiondir/$key.predict") ;
	my %result ;
	foreach my $line (@lines) {	
		if ($line =~ m/^(.*?)\s+(\d+)\s+(\d+)\s+(.)\d+/)  { 
			$result{$1}{start} 	= $2 ; 
			$result{$1}{end} 	= $3 ; 
			$result{$1}{strand} = $4 ; 
		} 
	}
	return %result ;
}

sub glimmer_2_GeneTable {
	# initialize the GeneTable
	my $key = shift ;
	my %glimmer = read_glimmer_predict($key); 
	foreach my $orf (sort keys %glimmer) {
		$GeneTable{$key}{$orf}{region_name} = $key ;
		$GeneTable{$key}{$orf}{region_size} = $AOI_sizes{$key} ;
		$GeneTable{$key}{$orf}{orf}	 		= $orf ;
		$GeneTable{$key}{$orf}{gene_name} 	= $orf ;
		$GeneTable{$key}{$orf}{gene_start} 	= $glimmer{$orf}{start} ;
		$GeneTable{$key}{$orf}{gene_end} 	= $glimmer{$orf}{end} ;
		$GeneTable{$key}{$orf}{gene_strand}	= $glimmer{$orf}{strand} ;
		$GeneTable{$key}{$orf}{gene_color} 	= $conf{defaultColor} ;
		$GeneTable{$key}{$orf}{function} 	= '' ;
		$GeneTable{$key}{$orf}{motifs} 		= '' ;
		$GeneTable{$key}{$orf}{annotation} 	= '' ;
		$GeneTable{$key}{$orf}{alignment_url} 	= '' ;
		$GeneTable{$key}{$orf}{alignment_filename} 	= '' ;
		$GeneTable{$key}{$orf}{dna}         = '';
		$GeneTable{$key}{$orf}{protein} 	= '';
		$GeneTable{$key}{$orf}{comment} 	= 'Glimmer ';
		$GeneTable{$key}{$orf}{blastjson} 	= '';
		# swap start end to make it uniform
		if ($GeneTable{$key}{$orf}{gene_start} > $GeneTable{$key}{$orf}{gene_end}) {
			my $tmp_start = $GeneTable{$key}{$orf}{gene_start} ;
			$GeneTable{$key}{$orf}{gene_start} = $GeneTable{$key}{$orf}{gene_end} ;
			$GeneTable{$key}{$orf}{gene_end} = $tmp_start ;
		}	
	}
}


sub GeneTable_2_ptt {
	#Location	Strand	Length	PID	Gene	Synonym	Code	COG	Product
	#1..1362	+	1362	-	dnaA	spr0001	-	-	
	my $outfilename = shift ;
	my @headers = ('region_name','region_size','gene_name','gene_start','gene_end','gene_strand','gene_color','function','motifs','annotation');
	foreach my $key (sort keys %GeneTable) {
		my @lines = join "\t", ('Location','Strand','Length','PID','Gene','Synonym','Code','COG	Product') ;
		foreach my $orf (sort keys %{$GeneTable{$key}}) {
			my $start = $GeneTable{$key}{$orf}{gene_start} ;
			my $end = $GeneTable{$key}{$orf}{gene_end} ;
			if ($start>$end) {my $tmp = $start; $start = $end ; $end=$tmp; }
			my @row = "$start..$end\t$GeneTable{$key}{$orf}{gene_strand}" ;
			push @row, $end-$start ;
			push @row, "-\t$orf\t$GeneTable{$key}{$orf}{gene_name}\t-\t-\t";
			push @lines, join "\t", @row ;
		}
		anne_files::write_lines("$sessiondir/$key.ptt", @lines);	
	}
}	


sub bacteriocin_hmmsearch {
	my $key = shift ;
	my $faa_file	= "$sessiondir/$key.faa" ; # the proteins in fasta format
	my $hmm_result_file 	= "$sessiondir/$key.bacteriocin_hmmsearch" ; # the Pfam search results
	unlink ($hmm_result_file) if (-e $hmm_result_file) ; # remove old result file because pfam_scan does not do this job
	my $command = "hmmsearch --tblout $hmm_result_file --domT 25 --cpu $cpu $program_dir/$conf{HMMs_folder}/bagel4_bacteriocin_hmms $faa_file >$sessiondir/null";
	system($command) ;
	print "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
	print "$command\n";
	print "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}


sub bacteriocin_hmm_2_GeneTable {
	# add hmmsearch results (4 columns) to GeneTable, and add the pfam/hmm description
	# orf00004	Barstar	PF01337	Barstar (barnase inhibitor)
	my ($key, $color, $inputfile) = @_ ;
	my @lines = anne_files::read_lines($inputfile) ;
	my @result ;
	foreach my $line (@lines) {
		if ($line !~ m/^\#/) {  # does not start with #
			my @items = split /\s+/, $line ;
			if (scalar @items > 4) {
				my ($pfam) = $items[3] =~ /(.*)\.\d+/ ;
				$hmm_descriptions{$pfam}{Description} = 'NA' if (!defined($hmm_descriptions{$pfam}{Description})) ;
				my $orf = $items[0] ;
				$GeneTable{$key}{$orf}{gene_name} 	= $items[0] ;
				$GeneTable{$key}{$orf}{gene_color} 	= $color ;
				$GeneTable{$key}{$orf}{function} 	.= $items[2].'; ' ;
				$GeneTable{$key}{$orf}{motifs} 		.= "<a href=http://pfam.xfam.org/family/$pfam>$pfam</a>".'; ' ;
				$GeneTable{$key}{$orf}{annotation} 	.= $hmm_descriptions{$pfam}{Description}.'; ' ;
				$GeneTable{$key}{$orf}{comment} 	.= "bacteriocin_hmm:$pfam ";
			}
		}
	}
}

sub blast_2_GeneTable {
	# add blast results to GeneTable, and add the description
	# line example: orf00022	Enterocin_X_chain_beta	53.19	47	22	0	5	51	1	47	5e-13	50.8
	my ($key, $inputfile, $db) = @_ ;
	print "=============> Webblast for $key db=$db file:$inputfile\n";
	my $inputhtml = $inputfile ;
		$inputhtml =~ s/$tmp_folder/$href_bagel4/ ;
	my @result ;
	my @lines = anne_files::read_lines($inputfile) ;
	my @json_filenames ; 
	foreach my $line (@lines) {
		if ($line !~ m/^\#/) {  # does not start with #
			my @items = split /\s+/, $line ;
			if (scalar @items > 4) {
				my $orf = $items[0] ;
				#my $json_file = "$key.$orf.blast_$db.json" ;
				my $json_file = "$key.$orf" ;
				$GeneTable{$key}{$orf}{gene_name} 	= $items[1] ;
				$GeneTable{$key}{$orf}{gene_color} 	= $conf{bacteriocinIColor} ;
				$GeneTable{$key}{$orf}{function} 	.= $items[1] ;
				$GeneTable{$key}{$orf}{annotation} 	.= "Evalue=$items[10] match=$items[2]%;" ;
				$GeneTable{$key}{$orf}{comment} 	.= "Bacteriocin blast";
				#$GeneTable{$key}{$orf}{blastjson} 	= $json_file;
				# make alignment HTML file and JSON for the GeneTablehtml report
				#if ($db eq '1' or $db eq '2') {  # add alignment file for bacteriocins of class I and II
				if ($db =~ m/(1|2)/) {  # add alignment file for bacteriocins of class I and II
					$GeneTable{$key}{$orf}{blastjson} 	= "$json_file.blast_$db.json";
					# add the webblast for detailed output of each blast hit orf
					anne_files::write_string("$sessiondir/$key.$orf.faa", ">$orf\n$GeneTable{$key}{$orf}{protein}\n");
					push @json_filenames, "$json_file.blast_$db.json";
					my $command = "$program_dir/bagel4_webblast.pl -s $sessiondir -i $key.$orf.faa -db '011' -o $json_file" ;  
					print "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
					print "$command\n";
					print "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n";
					system($command) ;					
					
					
				}	
			}
		}
	}
	anne_files::write_lines("$inputfile.filenames.table",@json_filenames) ;
	anne_files::print_lines(@json_filenames) ;
}


sub blast_bacteriocin_db {
	my $key = shift ;
	my $faa_file = "$sessiondir/$key.faa" ;
	my $blastresultfile = "$sessiondir/$key.bacteriocin.blast" ;
	my 	$tmp = "blastp -outfmt 6 -query $faa_file -db $program_dir/$conf{bacteriocinI_proteins}   -max_target_seqs 1 -num_threads $cpu -evalue $conf{blast_evalue_bacteriocinI} -out $blastresultfile".'_1' ;   system($tmp) ;
		print "=====> blastp db_1=$tmp\n" ;
		blast_2_GeneTable($key, $blastresultfile.'_1', '1');
	   
		$tmp = "blastp -outfmt 6 -query $faa_file -db $program_dir/$conf{bacteriocinII_proteins}  -max_target_seqs 1 -num_threads $cpu -evalue $conf{blast_evalue_bacteriocinII} -out $blastresultfile"."_2" ; 	system($tmp) ;
		blast_2_GeneTable($key, $blastresultfile.'_2', '2');

		$tmp = "blastp -outfmt 6 -query $faa_file -db $program_dir/$conf{bacteriocinIII_proteins} -max_target_seqs 1 -num_threads $cpu -evalue $conf{blast_evalue_bacteriocinIII} -out $blastresultfile"."_3" ; system($tmp) ;
		blast_2_GeneTable($key, $blastresultfile.'_3', '3');
}




sub translate_genes {
	my $key = shift ;
	my $genes_file = "$sessiondir/$key.genes";
	my $proteins_file =  "$sessiondir/$key.faa" ;
	my %genes = anne_genomics::fasta2hash($genes_file) ;
	my @lines ;
	foreach my $orf (sort keys %genes) {
		#print "$genes{$orf}\n";
		$GeneTable{$key}{$orf}{dna} = uc($genes{$orf}) ; 
		$GeneTable{$key}{$orf}{protein} = anne_genomics::translate($genes{$orf}) ; 
		$GeneTable{$key}{$orf}{protein} =~ s/\-//g ;
		push @lines, ">$orf\n".$GeneTable{$key}{$orf}{protein} ;	
	}
	anne_files::write_lines($proteins_file, @lines) ;
}

sub call_genes {
	# Orf calling on AOI using glimmer
	# Glimmer3.02
		# -o12  Set maximum overlap length to <n>.  Overlaps this short or shorter are ignored
		# -t12	Set threshold score for calling as gene to n
		# --orf_coords <filename>	Use <filename> to specify a list of orfs that should be scored separately, with no overlap rules
		# --separate_genes  <sequence-file> is a multifasta file of separate genes to     be scored separately, with no overlap rules
		# --no_indep predict more short genes
	my $AOI_key = shift ;
	my $AOI_fasta = shift ;
	make_glimmer_model($AOI_fasta) ;
	my $command = "$conf{glimmerpath}/glimmer3 -o 12 -g $conf{glimmer_min_gene_len} -t 12 --no_indep --linear $AOI_fasta $sessiondir/GLIMMER.icm $sessiondir/$AOI_key >>$sessiondir/glimmer.log";
	system($command) ;
	# Extract the genes
	$command = "$conf{glimmerpath}/extract -t $AOI_fasta $sessiondir/$AOI_key.predict > $sessiondir/$AOI_key.genes";
	system($command);	
}


sub make_glimmer_model {
	# Extract the training sequences from the genome file
	# Build the icm from the training sequences
	my $AOI_fasta = shift ;
	print "============== make glimmer model ================\n";
	my $command = "$conf{glimmerpath}/long-orfs -A atg -n -t 1.15 $AOI_fasta $sessiondir/GLIMMER.longorfs > /dev/null" ; #longorf only ATG allowed
	system($command) ; 
	$command = "$conf{glimmerpath}/extract -t $AOI_fasta $sessiondir/GLIMMER.longorfs > $sessiondir/GLIMMER.train";
	system($command);	
	$command = "$conf{glimmerpath}/build-icm -r $sessiondir/GLIMMER.icm < $sessiondir/GLIMMER.train";
	system($command);
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
    }
    die $usage if (!$fna_file) ;
}


