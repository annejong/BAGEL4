#!/usr/bin/env perl

#  Wrapper for BAGEL4 program 

# 	Anne de Jong
#	University of Groningen
#	the Netherlands
#	anne.de.jong@rug.nl
#
#  October 2017
#
#	This is the wrapper for the main pipeline of BAGEL4
#	this routine optional unzip the file and subsequently scans all the files and folders
#	From the webserver it just scans and analyse all the files, 
#	but as stand-alone the folder content can be filtered using a regular expression for the file names
#


use strict ;
use warnings ; 
use Data::Dumper ;
use File::Basename;
use lib "/data/molgentools/lib";
use anne_files ;
use lib "/data/bagel4/lib" ;
use bagel4_functions ;
use JSON;


my $program_dir = dirname($0) ;
my $sessiondir = "/data/bagel4/test" ;
my $query ;
my $webserver = 0 ;
my $rnaseq_example = 0 ;
my $rnaseq_bedgraph = '';
my $regex = '*';
my $usage = "./bagel4_wrapper.pl
				-s Sessiondir [default=$sessiondir]
				-query queryfolder with all the fasta files
				-webserver [default=0; Meaning that the analysis will de done on the folder from -i, else a query from the webserver is expected] 
				-r regular expression for selecting files from a folder [default= * :meaning all files ]
				-rnaseq_example [0|1] load the RNA-Seq example data
				-full full functional annotation of the genome context using PFAM 
		e.g.  /data/bagel4/bagel4_wrapper.pl -s /tmp/bagel4/test -query /var/genomes/g2d_mirror/Streptococcus_pneumoniae_D39 -r '*.fna'		
		e.g.  /data/bagel4/bagel4_wrapper.pl -s /tmp/bagel4/test -query /data/bagel4/examples -r 'pl*.fna'		
" ;

&parseparam() ;

my %conf = bagel4_functions::read_conf("$program_dir/bagel4.conf") ;
my $command = "$program_dir/bagel4_conf_2_json.pl -s $sessiondir" ;
system($command);

my $rnaseq_example_file = '/data/www/bagel4/php/examples/example.query.BedGraph';

# -----------------------------------------------------------------------------------  start main pipeline ---------------------------------------------------------------------------


mkdir $sessiondir ;

my %sessionprogress_new ; # new sessionprogress: global, query
$sessionprogress_new{global} = '';
$sessionprogress_new{query} = '';


# 1. Prepare the query folder
	add_sessionprogress('Checking input data', 'query') ;
	my @filenames ;		# old: all filenames
	if ($query =~ m/WGS\:(GCF_.*)\=/) { # Input is a WGS file and should be downloaded from NCBI
		@filenames = WGS_download($1) ;
	} elsif ($query =~ m/g2d_mirror.*\.fna$/) { # input is a single file (e.g. genome from the g2d mirror database selected from the webserver)
		@filenames = $query ; 
	} else { 
		print "=====================================>\n";	
		print "query=$query\n";
		print "regex=$regex\n";
		print "=====================================>\n";	

		my @all_filenames  = anne_files::get_files_from_subdirs_v2("$sessiondir/queryfolder", $regex) ; 
		foreach my $filename (@all_filenames) {
			if (uc($filename) =~ m/BEDGRAPH$/) { 
				 $rnaseq_bedgraph = $filename ;
			} else {
				push @filenames, $filename ;
			}
		}	
		if (check_php_code_injection(@filenames)) {
			&anne_files::write_log("$sessiondir/resulttable2.html","\n============ PHP code in upload is not allowed =============\n",'true');	
			&anne_files::write_log("$sessiondir/sessionstop","\n============ PHP code in upload is not allowed =============\n",'true');
			exit();
		}
	}	
	&anne_files::write_log("$sessiondir/BAGEL_wrapper.log","Scanning $query for files: $regex WEBSERVER=$webserver",'true');
	$sessionprogress_new{global} = 'Number of valid input files: '.(scalar @filenames).'<br>';

# 2. Scan the files and analyse them using the BAGEL pipeline
	my %filenames_db ;  # ALL filenames and entry names
	files_from_input(); # Fills %filenames_db. # scan all files and take the multiple fasta files, also size cutoff is done here
	my $InputfilesCount = scalar(keys %filenames_db) ;
	$sessionprogress_new{global} .= "Number of valid input DNA sequences: $InputfilesCount<br>";
	print "queryfolder=$query,  regex=$regex InputfilesCount=$InputfilesCount\n";	
	
	
	
	my $sequences_count = 0 ;
	my $DNAanalyzedSizeTotal =0 ;
	my $debug = 1 ;  # for debugging but not when run is from the webserver
	$debug = 0 if ($webserver) ;
	
	foreach my $key_db (sort keys %filenames_db) {
		my $queryname = $filenames_db{$key_db}{queryname} ;

		# some statistics
		$sequences_count++ ;
		$sessionprogress_new{query} = "Analyzing sequence <b>$queryname</b> [$sequences_count/$InputfilesCount] <br> " ;  # here we reset the web report for each queryname
		my $DNAanalyzedSize = length(anne_genomics::get_dna_from_fasta("$sessiondir/$queryname.fna")) ;
		$DNAanalyzedSizeTotal += $DNAanalyzedSize ;
		# old &anne_files::write_log("$sessiondir/BAGEL_wrapper.log","Analyzing nr $sequences_count/$InputfilesCount ==> $queryname",'true');
		
		
		# AOI identification on hmm's  tables/AOI_identification_hmm_rules
		add_sessionprogress('&emsp;AOI identification on the basis of motifs', 'query') ;
		translate_six_frames($queryname) ; # translate the query in the 6 possible reading frames
		system_command("$program_dir/bagel4_AOI_identification_hmm_rules.pl -s $sessiondir -dna $sessiondir/$queryname.fna -queryname $queryname") ;
		
		# AOI identification on blast
		add_sessionprogress('&emsp;AOI identification on the basis of blast', 'query') ;
		system_command("/data/bagel4/bagel4_AOI_identification_blast_bacteriocins.pl -s $sessiondir -translate true -dna $sessiondir/$queryname.fna -queryname $queryname") ;
	 
		# combine overlapping AOI's from hmmsearch and blast: AOI.blast.table AND AOI.hmm.table ==> $queryname.AOI.table
		system_command("/data/bagel4/bagel4_AOI_merge.pl -s $sessiondir -dna $sessiondir/$queryname.fna -queryname $queryname") ;
		

		# check if the AOI.table is empty
		if (-z "$sessiondir/$queryname.AOI.table") { 
			add_sessionprogress('&emsp;No AOIs found', 'query') ;
		} else {	
			add_sessionprogress("<font style=color:#0f9659>AOI(s) found in <b>$queryname</b></font>", 'global') ;
			
			# annotate the identified AOI's
			add_sessionprogress('&emsp;&emsp;Annotation of the identified AOIs', 'query') ;
			system_command("/data/bagel4/bagel4_AOI_annotation.pl -s $sessiondir -model false -dna $sessiondir/$queryname.AOI.fna -queryname $queryname")	;
			
			# make RNA-Seq tracks
			if (-e $rnaseq_bedgraph) { add_RNA_Seq_track($queryname, $rnaseq_bedgraph) ;} # RNA-Seq track is uploaded
			if ($rnaseq_example) { add_RNA_Seq_track($queryname, $rnaseq_example_file) ;} # RNA-Seq track from example data
			
			
			# promoter prediction to each AOI of the multi fasta
			add_sessionprogress('&emsp;&emsp;Promoter prediction', 'query') ;
			system_command("/data/bagel4/bagel4_promoterprediction.pl -s $sessiondir -dna $sessiondir/$queryname.AOI.fna") ;
			
			# terminator prediction to each AOI of the multi fasta
			add_sessionprogress('&emsp;&emsp;Transcription terminator prediction', 'query') ;
			system_command("/data/bagel4/bagel4_transterm.pl -s $sessiondir -dna $sessiondir/$queryname.AOI.fna") ;  # depends on .ptt file
			
			# JSON for Promoters and Terminators. make a json of the multifasta. Get the AOI's names from the fasta file headers
			add_sessionprogress('&emsp;&emsp;Create JSON files', 'query') ;
			system_command("/data/bagel4/bagel4_promoter_term_2_json.pl -s $sessiondir -queryname $queryname") ;

			# JSON for Genes. make a json of the multifasta. Get the AOI's names from the fasta file headers
			system_command("/data/bagel4/bagel4_GeneTable_2_json.pl -s $sessiondir -table $queryname.GeneTables") ;
		}
	}

# 3. Generate the overview table 	
	add_sessionprogress('Make overview table', 'global') ;
	make_OverviewGeneTables();

	print "Number of input files = $InputfilesCount\n";

# 4. Write the filenames to JSON
	anne_files::write_lines("$sessiondir/00.OverviewFilenames.json",(encode_json \%filenames_db)) ;

# tell the webserver that the session is finished
	&anne_files::write_log("$sessiondir/sessionstop","\n============ Analysis done =============\n",'true');

# Finally clean non useful or empty files
system_command("/data/bagel4/bagel4_clean_tmpfiles.pl -s $sessiondir") ;

# -----------------------------------------------------------------------------------  functions ---------------------------------------------------------------------------

sub system_command {
	my $command = shift ;
	print "+++++++++++++++++++++++++++++++++++++++++++\n$command \n+++++++++++++++++++++++++++++++++++++++++++\n";
	system($command) if (!$debug) ;
	anne_files::append_lines("$sessiondir/system_command.log", $command) ;
}

sub translate_six_frames {
	# translate DNA in all 6 frames and save to $queryname.six_frames.fna
	my $queryname = shift;
	my $DNA = anne_genomics::get_dna_from_fasta("$sessiondir/$queryname.fna") ;
	my %sixframes = bagel4_functions::six_frames($DNA) ;
	my @lines ;
	foreach my $frame (sort keys %sixframes) {
		push @lines, ">frame_$frame";
		push @lines, $sixframes{$frame} ;
	}	
	anne_files::write_lines("$sessiondir/$queryname.six_frames.fna", @lines) ;
}

sub WGS_download {
	my $assembly_genome = shift ;
	# read the summary file containing the ftp location
	my @lines = anne_files::read_lines("$program_dir/tools/00.latest_assembly_summary.txt");
	my $ftp ;
	my $genomename = '';
	foreach my $line (@lines) {
		if ($line =~ m/$assembly_genome/) {
			my @items = split /\t/, $line ;
			$ftp = $items[19] ;
			$genomename = ( split '/', $ftp )[ -1 ];
			system_command("wget $ftp/$genomename"."_genomic.fna.gz -O $sessiondir/$genomename.fna.gz") ;
			system_command("gunzip $sessiondir/$genomename.fna.gz") ;
			
			last;
		}
	}	
	return "$sessiondir/$genomename.fna" ;
}


sub add_sessionprogress {
	my ($str,$field) = @_ ;
	$sessionprogress_new{$field} .= $str.'<br>';
	#push @sessionprogress, $str.'<br>' ;
	anne_files::write_string("$sessiondir/sessionprogress",$sessionprogress_new{global}.'<hr>'.$sessionprogress_new{query}) ;
	
}



sub add_RNA_Seq_track {
	# 1. read the tracks and save for each queryname.AOI_XX.bedgraph.json
	my $queryname = shift ;
	my $bedgraph = shift ;
	my @lines = anne_files::read_lines("$sessiondir/$queryname.AOI.table");
	foreach my $line (@lines) {
		my @items = split /\t/, $line ;
		if (scalar @items>2) {
			my $command = "/data/bagel4/bagel4_bedgraph_2_json.pl -s $sessiondir -bed $bedgraph -queryname $items[0] -start $items[2] -end $items[3]" ;
			print "$command\n";
			system($command) if (!$debug) ;
		}	
	}
}



sub make_OverviewGeneTables {  # from GSEA-Pro
	my %OverviewTable ;
	
	my @html =  anne_files::read_lines("$program_dir/tables/GeneTable_script_header.html"); # js for hide/show detailed tables
	push @html, anne_files::read_lines("$program_dir/tables/GeneTable.css"); # css for detailed tables

	# get all the AOI tables of all DNA fragments
	my @table ;
	foreach my $key_db (sort keys %filenames_db) {
		my $queryname = $filenames_db{$key_db}{queryname} ;
		push @table, anne_files::read_lines("$sessiondir/$queryname.AOI.table") ;  # add the AOIs of each DNA
	}
	anne_files::print_lines(@table) ;
	my $AOI_count = scalar @table ;
	
	# ------------------ Run Summary table -------------------
	my $count_querynames = scalar(keys %filenames_db) ;
	my $count_filenames = scalar(@filenames) ;
	push @html,  "<TABLE id=SummaryTable width=50%>\n\t<tr><th>Run summary</th><th></th></tr>" ;
	push @html,  "\t<tr><td>Number of files analyzed</td><td>$count_filenames</td></tr>";
	push @html,  "\t<tr><td>Number of DNA fragments analyzed</td><td>$count_querynames</td></tr>";
	push @html,  "\t<tr><td>Total bases in all DNA</td><td>$DNAanalyzedSizeTotal</td></tr>";
	push @html,  "\t<tr><td>Number of AOI's <b>A</b>reas <b>O</b>f <b>I</b>nterest)</td><td>$AOI_count</td></tr>";
	push @html,  "</TABLE><br>\n";
	$OverviewTable{SummaryTable}{filecount} = $count_filenames ;
	$OverviewTable{SummaryTable}{querycount} = $count_querynames ;
	$OverviewTable{SummaryTable}{basecount} = $DNAanalyzedSizeTotal ;
	$OverviewTable{SummaryTable}{AOIcount} = $AOI_count ;
	
	# ------------------ URL + parameters to GeneTableGraphs.html ------------------
	my @sessionpath = split /\//, $sessiondir ;
	my $sessionID  = $sessionpath[-1];
	my $url_genetable_graph = "http://bagel4.molgenrug.nl/php/GeneTableSingleGraph.html?sessionID=$sessionID&AOI_count=$AOI_count&queryname=";

	# ------------------ AOI summary ------------------
	my %OverviewTableNew ;
	push @html, "<TABLE id=ResultsTable width=80%>\n\t<tr><th>AOI</th><th>start</th><th>end</th><th>Class</th><th>Fasta header</th></tr>";
	foreach my $line (@table) {
		$line =~ s/\|/; /g;  # replace the pipe by ;space
		my @items = split /\t/, $line ;
		if (scalar @items>4) {
			my $queryname = $items[0] ;
			
			my $querybasename = $queryname =~ s/\.AOI.*//gr ;
			my @key_filenames_db = grep { $filenames_db{$_}{queryname} eq $querybasename } keys %filenames_db;
			my @row = "<a href=$url_genetable_graph$queryname target=_blank><b>$items[0]</b></a>" ;
			push @row, ($items[2],$items[3],$items[5], $filenames_db{$key_filenames_db[0]}{fullname}) ;
			push @html, "\t<tr><td>".(join '</td><td>', @row)."</td></tr>\n" ;
			my %TableRecord ;
			$TableRecord{AOI} = "<a href=$url_genetable_graph$queryname target=_blank><b>$items[0]</b></a>" ;
			$TableRecord{start} = $items[2] ;
			$TableRecord{end} = $items[3] ;
			$TableRecord{class} = $items[5] ;
			$TableRecord{fullname} = $filenames_db{$key_filenames_db[0]}{fullname} ;
			$OverviewTableNew{$queryname}{AOI} = "<a href=$url_genetable_graph$queryname target=_blank><b>$items[0]</b></a>" ;
			$OverviewTableNew{$queryname}{start} = $items[2] ;
			$OverviewTableNew{$queryname}{end} = $items[3] ;
			$OverviewTableNew{$queryname}{class} = $items[5] ;
			$OverviewTableNew{$queryname}{fullname} = $filenames_db{$key_filenames_db[0]}{fullname} ;
			push @{$OverviewTable{ResultsTable}}, \%TableRecord ;
		}
	}
	push @html, '</TABLE>' ;
	anne_files::write_lines("$sessiondir/00.OverviewGeneTables.html",(@html)) ;
	anne_files::write_lines("$sessiondir/00.OverviewGeneTables.json",(JSON->new->pretty->utf8->encode(\%OverviewTable))) ;
	anne_files::write_lines("$sessiondir/00.OverviewGeneTablesNew.json",(JSON->new->pretty->utf8->encode(\%OverviewTableNew))) ;
}




sub OLD_Table_2_HTML {
	my $tab_table = shift ;
	my $html_table = shift ;
	my @lines = anne_files::read_lines($tab_table) ;
	my @html = anne_files::read_lines("$program_dir/table/jquery_sortable_table_header.html") ;
	push @html, "<TABLE id=myTable>";
	my $header = 1 ;
	foreach my $line (@lines) {
		my @elements = split /\t/, $line ;
		if ($header) { # header
			$header = 0 ;
			push @html, "<thead>\n\t<tr>" ;
			foreach (@elements) { push @html, "<th>$_</th>"; }
			push @html, "</tr>\n</thead><tbody>\n" ;
		} else { # body
			push @html, "\t<tr>" ;
			foreach (@elements) { push @html, "<td>$_</td>"; }
			push @html, "</tr>\n" ;
		}	
	}	
	push @html, "</tbody></TABLE>\n" ;
	anne_files::write_lines($html_table, @html) ;
}


sub fasta2hash {
	# the header ">.." is the key
	# see also g2dfasta2hash
	my $file = shift ;
	my %result ;
	my $seq = '';
	my $key = '' ;
	my @lines = anne_files::read_lines($file) ;
	foreach my $line (@lines) {
		if ($line =~ m/^\>(.*)/) { # new record, add the last to the hash
			$result{$key} = $seq if ($key ne '') ;
			$seq = ''; 
			$key = $1 ;
			#$key = anne_misc::trim($1) ; # remove leading and taling space
			#if ($key =~ /(.*?)\s/) { $key = $1; }
		} else {
			$seq .= $line ;
		}
	}
	# don't forget the last record
	$result{$key} = $seq if ($key ne '') ;
	return %result ;
}

sub files_from_input {
	# input can be multiple multiple-fasta files 
	my $key_db = 0 ;
	my $folders = "/data/g2d_mirror/|$sessiondir/queryfolder/" ;
	my @lines = "ID\tfilename\tqueryname\tlength" ;
	my $ID = 0 ;
	foreach my $file (@filenames) {
		# check on multifasta
		#my %fasta = anne_genomics::fasta2hash($file);
		my %fasta = fasta2hash($file);
		foreach my $key (keys %fasta) {
			if (length($fasta{$key}) > $conf{minimum_region}) {
				my $queryname = $key;
				if ($queryname =~ m/(.*?)\s/) { $queryname = $1; }
				if ($queryname =~ m/^(.*)( |)/g) { $queryname = $1; } # take only the info before the |
				$queryname = &anne_files::properfilename($queryname).".$key_db" ;  # remove unwanted chars AND add unique number to the queryname to enable the user to use same names in different files
				anne_files::write_string("$sessiondir/$queryname.fna", ">$queryname\n$fasta{$key}\n") ;
				$key_db++;
				$file =~ s/$folders//;
				
				$filenames_db{$key_db}{filename} = $file ;
				$filenames_db{$key_db}{fullname} = $key ;
				$filenames_db{$key_db}{queryname} = $queryname ;
				$filenames_db{$key_db}{length} = length($fasta{$key}) ;
				$ID++;
				push @lines, "$ID\t$filenames_db{$key_db}{filename}\t$filenames_db{$key_db}{queryname}\t$filenames_db{$key_db}{length}";
			}	
		}
	}
	anne_files::write_lines("$sessiondir/00.all_contigs.table", @lines);
	anne_files::write_lines("$sessiondir/00.filenames_db.json",(JSON->new->pretty->utf8->encode(\%filenames_db))) ;
}


sub check_php_code_injection {
	# will return true is PHP code is found one of the files
	my @filenames = @_ ;
								   
	my $result = 0 ;
																		
	foreach my $file (@filenames) {
		print "Checking for PHP code: $file\n";
		open(FILE,$file);
		my $content = <FILE>;
		if ($content =~ /\<\?php/) {  # check if the file contains php code
			$result = 1 ;
			last ;
		}	
		close(FILE);
	}
	return $result ;
}


sub parseparam {
    my $var ;
    my @arg = @ARGV ;
    while(@arg) {
        $var = shift(@arg) ;
		die $usage if ($var eq '-h' or $var eq '--help') ;
		$sessiondir		= shift(@arg) if($var eq '-s') ;
		$query			= shift(@arg) if($var eq '-query') ;
		$webserver		= shift(@arg) if($var eq '-webserver') ;
        $rnaseq_example	= shift(@arg) if($var eq '-rnaseq_example') ;
        $regex			= shift(@arg) if($var eq '-r') ;
    }
    die $usage if (!$query) ;
	$sessiondir =~ s/\/$// ; #remove the last / from sessiondir to make it universal
}

