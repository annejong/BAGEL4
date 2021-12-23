#!/usr/bin/env perl

use strict ;
use warnings ;
use lib "/data/bagel4/lib" ;
use bagel4_functions ;
use bagel4_blast ;
use lib "/data/molgentools/lib";
use anne_files ;
use File::Basename;
use LWP::Simple;


my $sessiondir = '/data/bagel4/test';
my $query ;
my $cpu = 4 ;
my $databases = '011' ;
my $max_targets = 1 ;
#my $evalue = 0.00001 ;
my $evalue = 0.1 ;
my $program_dir = dirname($0) ;
my $AOIoutput = '' ;

my $usage = "options:
				-s Sessiondir [default=current folder]
				-i query
				-db 011  bacteriocinI = +1 ; bacteriocinII = +10 ; bacteriocinIII = +100 ; 
				-m max_targets [default=1]
				-o optional for call from AOI annotation
		e.g.  ./bagel4_webblast.pl -s /data/bagel4/test -i webblast_query.faa 
\n - Anne de Jong - \n" ;



&parseparam() ;

my %conf = bagel4_functions::read_conf("$program_dir/bagel4.conf") ;

# ---------------------------------------------------------- main -------------------------------------------------------------------------------

# 1. Read the protein databases
	my %proteins_db = (	anne_files::read_table_to_hash("$program_dir/db_proteins/bagel4_class1bacteriocin_db.txt"),
						anne_files::read_table_to_hash("$program_dir/db_proteins/bagel4_class2bacteriocin_db.txt"),
						anne_files::read_table_to_hash("$program_dir/db_proteins/bagel4_class3bacteriocin_db.txt") ) ;

	
			my $tmp = "/data/molgentools/tools/dos2linux.pl $sessiondir/$query";
			print "$tmp\n";
			system($tmp) ;	
	my $Query_seq = getquerydna("$sessiondir/$query");
	


	
# 1. perform the blast
	my $db = "$program_dir/$conf{bacteriocinI_proteins}" ;
	for (my $i = 1; $i <=3; $i++) {
		if (substr($databases, 3-$i, 1) eq '1') {
			my $db ;
			$db = "$program_dir/$conf{bacteriocinI_proteins}" if ($i == 1) ;
			$db = "$program_dir/$conf{bacteriocinII_proteins}" if ($i == 2) ;
			$db = "$program_dir/$conf{bacteriocinIII_proteins}" if ($i == 3) ;
			my $tmp = "blastp -outfmt 5 -query $sessiondir/$query -db $db -max_target_seqs 5 -num_threads $cpu -evalue $evalue -out $sessiondir/my_blastresuts_".$i.".xml";
			print "$tmp\n";
			system($tmp) ;
			
			parse_outfmt_5("$sessiondir/my_blastresuts_".$i.".xml", "$sessiondir/webblast.bacteriocin_".$i.".html", "$AOIoutput.blast_$i.json") ;
		}	
	}
	

# 2. report if no result is found
if (!glob("$sessiondir/webblast.bacteriocin_*.html")) { 
	anne_files::write_string("$sessiondir/webblast.noresults.html", '<font color=blue><h2>No hits found</h2></font><br>') ; 
}	

# 3. tel webserver that run is finished
anne_files::write_string("$sessiondir/sessionstop", 'Session ready') if ($AOIoutput eq '') ;  # added by Anne: do not finish session if call is from AOI annotation


# ---------------------------------------------------------- functions -------------------------------------------------------------------------------

 sub getquerydna {
	my $infile = shift ;
	my @lines =  anne_files::read_lines($infile) ;
	my $result ;
	foreach my $line (@lines) {
		if ($line !~ m/^>/) { $result .= $line; }
	}
	$result =~ s/(\s|\x0d|\d|\|)//g;
	return $result ;
}	

 sub parse_outfmt_5 { # outputformat as XML
	my ($infile, $outfile, $JSONoutput) = @_ ;
	my @lines = anne_files::read_lines($infile) ;
	my $record = 0 ;
	my %blast ;
	my $key = 0 ;
	foreach my $line (@lines) {
		my $qlen=$1	if ($line =~ m/<BlastOutput_query-len>(.*)<\/BlastOutput/ ) ;
		$key++ if ($line =~ m/<Hit_id>(.*)<\/Hit/ ) ;
		$blast{$key}{Hit_id} = $1 			if ($line =~ m/<Hit_id>(.*)<\/Hit/ ) ;
		$blast{$key}{Hit_def} = $1 			if ($line =~ m/<Hit_def>(.*)<\/Hit/ ) ;
		$blast{$key}{Hsp_bit_score} = $1 	if ($line =~ m/<Hsp_bit-score>(.*)<\/Hsp/ ) ;
		$blast{$key}{Hit_accession} = $1 	if ($line =~ m/<Hit_accession>(.*)<\/Hsp/ ) ;
		$blast{$key}{Hsp_evalue} = $1 		if ($line =~ m/<Hsp_evalue>(.*)<\/Hsp/ ) ;
		$blast{$key}{Hsp_qseq} = $1 		if ($line =~ m/<Hsp_qseq>(.*)<\/Hsp/ ) ;
		$blast{$key}{Hsp_hseq} = $1 		if ($line =~ m/<Hsp_hseq>(.*)<\/Hsp/ ) ;
		$blast{$key}{Hsp_hitfrom} = $1 		if ($line =~ m/<Hsp_hit-from>(.*)<\/Hsp/ ) ;
		$blast{$key}{Hsp_hitto} = $1 		if ($line =~ m/<Hsp_hit-to>(.*)<\/Hsp/ ) ;
		$blast{$key}{Hit_len} = $1 			if ($line =~ m/<Hit_len>(.*)<\/Hit/ ) ;
		$blast{$key}{Hsp_midline} = $1 		if ($line =~ m/<Hsp_midline>(.*)<\/Hsp/ ) ;
		$blast{$key}{Hsp_query_from} = $1 	if ($line =~ m/<Hsp_query-from>(.*)<\/Hsp/ ) ;
		$blast{$key}{Hsp_query_to} = $1 	if ($line =~ m/<Hsp_query-to>(.*)<\/Hsp/ ) ;
		$blast{$key}{qlen} = $qlen ;
	}
	my @html ;
	my @json = '{"blasthits": [' ; 
	foreach my $key (sort keys %blast) {
		if (defined($blast{$key}{Hit_def})) {
			$blast{$key}{Hsp_midline} =~ s/\s/&nbsp;/g ;   # space to html-space
			$blast{$key}{Hit_accession} = 'NA' if (!defined($blast{$key}{Hit_accession})) ;
			my $firsthit = 1 ;
			my $leaderhighlighted = 0 ;
			if ($blast{$key}{Hit_id} =~ m/(.*);(.*)/) { 
				push @json, '{ "blasthit": "'.$key.'",' ;
				my $db_ID = $1 ; 
				my $hit_ID = $2 ;
				my $reference = '';  $reference = "<a href=$proteins_db{$db_ID}{'http://dx.doi.org/'} target=_blank>Reference</a>" if ($proteins_db{$db_ID}{'http://dx.doi.org/'} ne '') ;
				my %features ; 
				if ($firsthit) { 
					$firsthit = 0 ; 
					%features = bagel4_blast::get_uniprot_features("$program_dir/db_uniref90_xml/$proteins_db{$db_ID}{Uniprot_ID}.xml") ; 
				}
				
				
				
				
				
				
				
				$blast{$key}{Hsp_hseq}= substr($proteins_db{$db_ID}{Sequence}, 0, ($blast{$key}{Hsp_hitfrom}-1)).($blast{$key}{Hsp_hseq}).substr($proteins_db{$db_ID}{Sequence}, $blast{$key}{Hsp_hitto}, (length ($proteins_db{$db_ID}{Sequence})-$blast{$key}{Hsp_hitto}));
				$blast{$key}{Hsp_midline}= '&nbsp' x ($blast{$key}{Hsp_hitfrom}-1).($blast{$key}{Hsp_midline});
				# DEDUG by Anne: This adds all AA to the query when blasting multiple proteins  $blast{$key}{Hsp_qseq}= '&nbsp' x ($blast{$key}{Hsp_hitfrom}-$blast{$key}{Hsp_query_from}).substr($Query_seq, 0, ($blast{$key}{Hsp_query_from}-1)).($blast{$key}{Hsp_qseq}).substr($Query_seq, $blast{$key}{Hsp_query_to}, (length ($Query_seq)-$blast{$key}{Hsp_query_to}));
				$blast{$key}{Hsp_qseq}= '&nbsp' x ($blast{$key}{Hsp_hitfrom}-$blast{$key}{Hsp_query_from}).substr($Query_seq, 0, ($blast{$key}{Hsp_query_from}-1)).($blast{$key}{Hsp_qseq});
				
				
				
				# Auke generate a line with all the modified residues highlighted
				my @modifiedresidues;
				foreach my $feature_key (keys %features) { 
					if ($features{$feature_key}{type} eq 'modified residue') {
					push @modifiedresidues, $features{$feature_key}{position};
					}
				}
				my @sortedmodifiedresidues = sort { $b <=> $a } @modifiedresidues;  
				my $modifiedresidues = $blast{$key}{Hsp_hseq};
				foreach my $mr (@sortedmodifiedresidues){ 
				$modifiedresidues = substr($modifiedresidues, 0,($mr-1)).'<font style=color:red>*</font><font style=color:white#EEEEEE;>'.substr($modifiedresidues, $mr);
				}
				my $modifiedresidues2 = '<font style=color:#EEEEEE;>'.$modifiedresidues.'</font>';
				
				#being smart about gaps in the hseq sequence
				
				my @hseqgaps;
				my $offsetx = 0;
				my $blankspace2=0;
				my $resultx = index ($blast{$key}{Hsp_hseq},"-", $offsetx);
				while ($resultx != -1) {
					
					if ($blast{$key}{Hsp_hitfrom} < $blast{$key}{Hsp_query_from}){
						$blankspace2=($blast{$key}{Hsp_query_from}-$blast{$key}{Hsp_hitfrom});
						}
					push @hseqgaps, ($resultx+$blankspace2);
					$offsetx = $resultx + 1;
					$resultx = index($blast{$key}{Hsp_hseq}, "-", $offsetx);
					}
				
								
				
				# Auke solve the rare issue that query from is larger then hit from
				my $blankspace =0;
				
				if ($blast{$key}{Hsp_hitfrom} < $blast{$key}{Hsp_query_from}) {
						$blast{$key}{Hsp_midline}= '&nbsp' x ($blast{$key}{Hsp_query_from}-$blast{$key}{Hsp_hitfrom}).($blast{$key}{Hsp_midline});
						$blast{$key}{Hsp_hseq}= '&nbsp' x ($blast{$key}{Hsp_query_from}-$blast{$key}{Hsp_hitfrom}).$blast{$key}{Hsp_hseq};
						$blankspace = ($blast{$key}{Hsp_query_from}-$blast{$key}{Hsp_hitfrom});
						$modifiedresidues2 = '&nbsp' x ($blast{$key}{Hsp_query_from}-$blast{$key}{Hsp_hitfrom}).$modifiedresidues2;
						}
				
				
				
				#Auke make a line indicating bridges
				my %bridges;
				my $id = 0;
				foreach my $feature_key (keys %features) { 
					if ($features{$feature_key}{type} =~ m/(cross-link|disulfide bond)/) {
						$id++;
						$bridges{$id}{beginpos}=$features{$feature_key}{beginpos}+$blankspace;									
						$bridges{$id}{endpos}=$features{$feature_key}{endpos}+$blankspace;
					}
					
				}
				
				#Auke generating single bridging line
								
				my @bridgelines;
				my $bridgeline = ("1"x length($blast{$key}{Hsp_hseq}));
				my $brdigelinetest = $bridgeline;
				my $currentbridgeline ="";
				#push @html, $bridgeline;
				foreach my $id (sort {$bridges{$a}{beginpos} <=> $bridges{$b}{beginpos}} keys %bridges) {
					
					if ((substr $bridgeline, $bridges{$id}{beginpos}-1 , 1) eq 3) {  #if overlapping bridges save line and 
						push @bridgelines, $currentbridgeline;
						$bridgeline = ("1"x length($blast{$key}{Hsp_hseq})); #resetting the bridgeline
						}
					 
					$bridgeline = (substr $bridgeline, 0,$bridges{$id}{beginpos}-1)."2"."3"x ($bridges{$id}{endpos}-$bridges{$id}{beginpos}-1)."4".(substr $bridgeline, $bridges{$id}{endpos}, length($blast{$key}{Hsp_hseq})-$bridges{$id}{endpos});					
					$currentbridgeline=$bridgeline;
					}
				
				if ($bridgeline ne $brdigelinetest){ push @bridgelines, $bridgeline; }#prevent an empty bridgeline
				
				foreach my $iddd (@bridgelines){ #inserting the gap in the blast hit into the bridgeline
					foreach my $gapss (@hseqgaps){
						$iddd = (substr $iddd, 0, $gapss)."1".(substr $iddd, $gapss )
					}
				}
				foreach (@bridgelines) {
					s/2/└/g;
					s/3/─/g;
					s/4/┘/g;
					s/[1]/&nbsp;/g;
				}
				
				# Highlighting the leader /  signal peptide from UniProt data
				foreach my $feature_key (keys %features) { 
					if (($features{$feature_key}{type} eq 'propeptide') or($features{$feature_key}{type} eq 'signal peptide')) { $leaderhighlighted = 1 ;
						$blast{$key}{Hsp_hseq}  =~ s/nbsp//g;	
						$blast{$key}{Hsp_hseq} = (substr $blast{$key}{Hsp_hseq},0,$blankspace).'<font style=background-color:lightgrey;color:black>'.(substr $blast{$key}{Hsp_hseq},$blankspace,$features{$feature_key}{endpos}).'</font>'.(substr $blast{$key}{Hsp_hseq},($features{$feature_key}{endpos}+$blankspace)) ;
						$blast{$key}{Hsp_hseq} =~ s/\&/\&nbsp/g;
					}
				}	
				
				# Highlighting the leader from the database information but only if not done so before
				if ($proteins_db{$db_ID}{Leaderlength} && $leaderhighlighted ne 1 ) {
					$blast{$key}{Hsp_hseq}  =~ s/nbsp//g;
					$blast{$key}{Hsp_hseq} = (substr $blast{$key}{Hsp_hseq},0,$blankspace).'<font style=background-color:lightgrey;color:black>'.(substr $blast{$key}{Hsp_hseq},$blankspace,$proteins_db{$db_ID}{Leaderlength}).'</font>'.(substr $blast{$key}{Hsp_hseq},($proteins_db{$db_ID}{Leaderlength})+$blankspace) ;
					$blast{$key}{Hsp_hseq} =~ s/\&/\&nbsp/g;
				}

				# Highlighting Cysteins in lanthipeptides and sactipeptides
				if(substr($proteins_db{$db_ID}{Subclass},0,4) eq 'lant' or 'Sact'){$blast{$key}{Hsp_qseq} =~ s/C/\<font color=red\>C\<\/font\>/g ;}
				if(substr($proteins_db{$db_ID}{Subclass},0,4) eq 'lant' or 'Sact'){$blast{$key}{Hsp_midline} =~ s/C/\<font color=red\>C\<\/font\>/g ;}
				if(substr($proteins_db{$db_ID}{Subclass},0,4) eq 'lant' or 'Sact'){$blast{$key}{Hsp_hseq} =~ s/C/\<font color=red\>C\<\/font\>/g ;}
				
				
				#<Iteration_query-def> could be used to add input name to alignment
				# Make the html table for current record [anne]
				my @html_table;
				push @html_table, '<TABLE class=webblast>';
				push @html_table, "<tr><th></th><th align=center>Database hit with <b>$hit_ID</b> (Bit Score=$blast{$key}{Hsp_bit_score})</th></tr>";
				push @html_table, "<tr><td><br/></td><td></td></tr>";
				push @html_table, "<tr><td class=align_right><b>Query</b></td><td><font face=courier>$blast{$key}{Hsp_qseq}</font></td></tr>";
				push @html_table, "<tr><td></td><td><font face=courier>$blast{$key}{Hsp_midline}</font></td></tr>";
				push @html_table, "<tr><td class=align_right><b>$hit_ID</b></td><td><font face=courier>$blast{$key}{Hsp_hseq}</font></td></tr>";
				if (@bridgelines){ 
					foreach my $varrr (@bridgelines){ push @html_table, "<tr><td class=align_right>Bridges</td><td><font face=courier>$varrr</font></td></tr>" ; }
				}
				
				if (@modifiedresidues){ push @html_table, "<tr><td class=align_right>Modifications</td><td><font face=courier>$modifiedresidues2</font></td></tr>"; }
				push @html_table, "<tr><td><br/></td><td><font face=courier></font></td></tr>";
				if($proteins_db{$db_ID}{Subclass}){				push @html_table, "<tr><td class=align_right>Subclass</td><td>$proteins_db{$db_ID}{Subclass}</td></tr>";}
				if($proteins_db{$db_ID}{ORGANISM}){				push @html_table, "<tr><td class=align_right>Organism</td><td>$proteins_db{$db_ID}{ORGANISM}</td></tr>";}
				if($proteins_db{$db_ID}{'http://dx.doi.org/'}){	push @html_table, "<tr><td class=align_right>Literature</td><td>$reference</td></tr>";}
				if($proteins_db{$db_ID}{NCBI_ID}){				push @html_table, "<tr><td class=align_right>NCBI</td><td><a href='https://www.ncbi.nlm.nih.gov/protein/$proteins_db{$db_ID}{NCBI_ID}' target=_blank>$proteins_db{$db_ID}{NCBI_ID}</a></td></tr>";}
				if($proteins_db{$db_ID}{Uniprot_ID}){			push @html_table, "<tr><td class=align_right>UniProt</td><td><a href='http://www.uniprot.org/uniprot/$proteins_db{$db_ID}{Uniprot_ID}' target=_blank>$proteins_db{$db_ID}{Uniprot_ID}</a></td></tr>";}
				
				foreach my $feature_key (sort keys %features) { 
					push @html_table, "<tr><td></td><td>$features{$feature_key}{type} $features{$feature_key}{description} $features{$feature_key}{position}"; 
					if ($features{$feature_key}{endpos}){ push @html_table, "$features{$feature_key}{beginpos} - $features{$feature_key}{endpos} "; }
					push @html_table, "</td></tr>";
				}
				push @html_table,	"</TABLE>";
				push @html, @html_table ; # add the new hit to the main html
				my $html_record = join "", @html_table ;
				push @json, ' "htmlrecord": "'.$html_record.'"}';
				push @json, ',';
				

				
				
			}	
		}
	}
	
	pop @json; # remove the last comma
	push @json, ']}' ;
	anne_files::write_lines($outfile,@html) if (scalar @html > 3) ; 
	anne_files::write_lines("$sessiondir/$JSONoutput",@json) ; # added by Anne: call from the AOI annotation
}



sub parseparam {
    my $var ;
    my @arg = @ARGV ;
    while(@arg) {
        $var = shift(@arg) ;
        die $usage if ($var eq '-h' or $var eq '--help') ;
        $sessiondir 	= shift(@arg) if($var eq '-s') ;
		$query			= shift(@arg) if($var eq '-i') ;
		$databases			= shift(@arg) if($var eq '-db') ;
		$max_targets	= shift(@arg) if($var eq '-m') ;
		$evalue	= shift(@arg) if($var eq '-evalue') ;
		$AOIoutput	= shift(@arg) if($var eq '-o') ;
    }
    die $usage if (!$query) ;
}


