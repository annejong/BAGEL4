package anne_genomics ;

##############################################################
##      
## 	Anne de Jong
##	University of Groningen
##	the Netherlands
##	anne.de.jong@rug.nl
##
##############################################################
##
##   
##  2013 September, Copy from SERVER to NGS and removed all anne_files duplicated
##
##

use strict ;
use warnings ;
use Bio::Seq ;
use Bio::Tools::pICalculator;
use Bio::Tools::SeqPattern ;
use Bio::SearchIO ;
use Bio::SeqIO ;
use Bio::SeqIO::fasta;
use DBI ;
use lib "/usr/molgentools/lib";
use anne_files;
use anne_misc;

BEGIN {
    use Exporter ();
    @genomics::ISA         = qw(Exporter);
    @genomics::EXPORT      = qw();
    @genomics::EXPORT_OK   = qw(%genome %genomename %codontable);
}
use vars qw( %genome %genomename %codontable) ;

my %codontable ;
my $codontablefile ='/data/molgentools/genomics/codon_table.txt';
my $sessiondir = '.';

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------


# sub read_BAM {
	# my @sam ;
	# open BAM,"samtools view <bam_file> |";
	# while(<BAM>){
		# next if(/^(\@)/);  ## skipping the header lines (if you used -h in the samools command)
		# s/\n//;  s/\r//;  ## removing new line
		# @sam = split(/\t/);  ## splitting SAM line into array
	# }
	# chomp(@sam) ;
	# return @sam
# }

sub multi_2_singlefasta {
	# convert / split one multi entry FASTA file to single files and save them in sessiondir
	my $filename = shift ;
	my $sessiondir = shift ;
	my @lines = read_lines($filename);
	my $fastafile ;
	my @fasta ;
	my $firstline = 1 ;
	foreach my $line (@lines) {
		if ($line =~ m/^\>(.*)/) {  # Header found
			if ($firstline) {
				$firstline = 0 ;
			} else {	
				anne_files::write_lines("$sessiondir/$fastafile", @fasta) ;
			}
			$fastafile = $1 ;
			@fasta = ">$fastafile";
		} else {
			push @fasta, $line ;
		}
	}
	anne_files::write_lines("$sessiondir/$fastafile", @fasta) ; # write the last record
}



sub read_terminators {
	# Parse results of transtermHP
	my %result ;
	my @lines = anne_files::read_lines(@_) ;
	foreach my $line (@lines) {
		$line =~ s/(\n|\r|\x0d)//g;		 # remove line breaks etc
		my @items = split(/\s+/, $line) ;
		if (defined($items[3]) and $line !~ m/NONE/ and $items[2] =~ m/\d+/) {
			my $key = $items[1] ;
			$result{$key}{start} 	= $items[2] ; 
			$result{$key}{end} 		= $items[4] ; 
			$result{$key}{position}	= $items[2]+20 ; 
			$result{$key}{strand} 	= $items[5] ; 
			#$result{$key}{name} 	= $items[1] ; 
			$result{$key}{score} 	= $items[6] ; 
		}
	}
	return %result ;
}

sub gff_2_hash {
	my $gff_filename = shift ;
	my @lines = anne_files::read_lines($gff_filename);
	my %table ;
	my $key = 0 ;
	my $previous_locus = '';
	my @IDfields = ('Name','gene','locus_tag','GeneID','Note','product','protein_id','Dbxref') ;  # the fields to use from the ID column
	foreach my $line (@lines) {
		
		if ($line !~ m/^\#/) {  # line does not start with a #
			my @fields = split /\t/, $line ;
			$key++ ;  # use numbering as an unique identifier
			$table{$key}{NODE} 		= $fields[0] ; 	
			$table{$key}{genome} 	= $fields[0] ; 	
			if ($table{$key}{genome} =~ m/(.+)\..*/ ) { $table{$key}{genome} = $1 ; }  # remove version number
			$table{$key}{genome} 	=~ s/\s+//g ;  # remove spaces
			$table{$key}{type} 		= $fields[2] ; 	
			$table{$key}{start} 	= $fields[3] ; 	
			$table{$key}{end} 		= $fields[4] ; 	
			$table{$key}{strand}	= $fields[6] ;
			$table{$key}{length} 	= $table{$key}{end} - $table{$key}{start} +1 ;
			$table{$key}{RNA}		= '' ;
			my @IDs = split /\;/, $fields[8] ;		# the ID field
			foreach my $IDfield (@IDfields) { $table{$key}{$IDfield} = ''; } # be sure that the cells are defined 
			foreach my $ID (@IDs) {
				foreach my $IDfield (@IDfields) {
					if ( $ID =~ m/.*$IDfield=(.*)/ ) 	{ 
						my $clean_field = $1 ;
						$clean_field =~ s/(\:|\=)//g ;
						$table{$key}{$IDfield} = $clean_field ; 
						if ($table{$key}{type} eq 'gene' and $IDfield eq 'gene') { $table{$key}{gene} = $clean_field; } 
						elsif ($table{$key}{type} eq 'CDS'  and $IDfield eq 'protein_id') { $table{$key}{protein_id} = $clean_field; } 
						elsif ($table{$key}{type} =~ m/RNA/  and $IDfield eq 'product') { $table{$key}{RNA} = $clean_field; }
					}	
				}
			}
			if ($table{$key}{type} ne 'gene') {
				$table{$key}{locus_tag} = $previous_locus ;
			} else {
				$previous_locus = $table{$key}{locus_tag} ;
			}
			#print "$key;$table{$key}{genome};$table{$key}{gene};$table{$key}{locus_tag};$table{$key}{start};$table{$key}{end};$table{$key}{strand}\n";
		}
	}

	return %table ;
}

sub gff_remove_CDS {
	my %gff = @_ ;
	foreach my $key ( sort keys %gff) {	if ($gff{$key}{type} eq 'CDS' or $gff{$key}{type} eq 'region') { delete $gff{$key} ; } }
	return %gff ;
}

sub gff_2_hash_v3_old {
	# Parse the GFF annotation file to a hash
	# Can extract GFF 'type' such as RNA, gene, CDS, tRNA, rRNA; use regular expression e.g. gene|RNA 	
	# added: aug 2014
	my $gff_filename = shift ;
	my $gff_type = shift ; 
	#print "Reading gff file $gff_filename\n";
	my @lines = anne_files::read_lines($gff_filename);
	my %result ;
	foreach my $line (@lines) {
		if ($line !~ m/^\#/) {		# line does not start with #
			my @items = split /\t/ , $line ;
			if (@items>2) {
				if ( $items[2] =~ /($gff_type)/ ) {  # only take this type which match the $gff_type
					my $genome		= $items[0] ;
					my $type		= $items[2] ;
					my $ID	 = '' ;
					if ($items[8] =~ m/.*ID=(.*);Name=(.+?);/) { $ID = $1 ; }		# gene
					if ($items[8] =~ m/.*ID=(.*);Parent=(.+?);/) { $ID = $2 ; }	# protein
					if ( $ID ne '' ) {
						$result{$genome}{$type}{$ID}{Name} 		= $2 ;
						$result{$genome}{$type}{$ID}{locus_tag}	= $ID ;
						$result{$genome}{$type}{$ID}{genome} 	= $items[0] ;
						$result{$genome}{$type}{$ID}{contig} 	= $items[0] ;
						$result{$genome}{$type}{$ID}{start} 	= $items[3] ;
						$result{$genome}{$type}{$ID}{end} 		= $items[4] ;
						$result{$genome}{$type}{$ID}{strand} 	= $items[6] ;
						$result{$genome}{$type}{$ID}{length} 	= $result{$genome}{$type}{$ID}{end} - $result{$genome}{$type}{$ID}{start} +1 ;
						$result{$genome}{$type}{$ID}{gene}		= '';
						$result{$genome}{$type}{$ID}{old_locus_tag}	= '';
						$result{$genome}{$type}{$ID}{Parent} 	= '' ;
						$result{$genome}{$type}{$ID}{product}	= '';
						$result{$genome}{$type}{$ID}{note}		= '';
						$result{$genome}{$type}{$ID}{pseudo}	= 'false' ;
						my @descriptions = split ';', $items[8] ;
						$result{$genome}{$type}{$ID}{gene}		= $1 if ($descriptions[1]=~ m/Name=(.*)/ ) ;
						foreach my $description (@descriptions) {
							$result{$genome}{$type}{$ID}{locus_tag}	= $1 if ( $description =~ m/^locus_tag=(.*)/ ) ;
							$result{$genome}{$type}{$ID}{old_locus_tag}	= $1 if ( $description =~ m/^old_locus_tag=(.*)/ ) ;
							$result{$genome}{$type}{$ID}{gene} 		= $1 if ( $description =~ m/^gene=(.*)/ ) ;
							$result{$genome}{$type}{$ID}{Parent} 	= $1 if ( $description =~ m/^Parent=(.*)/ ) ;
							$result{$genome}{$type}{$ID}{product}	= $1 if ( $description =~ m/^product=(.*)/ ) ;
							$result{$genome}{$type}{$ID}{note}		= $1 if ( $description =~ m/^Note=(.*)/ ) ;
							$result{$genome}{$type}{$ID}{pseudo}	= $1 if ( $description =~ m/^pseudo=(.*)/ ) ;
						}	
					}	
				}	
			}	
		}		
	}
	return %result ;
}


sub gff_2_hash_v3 {
	# Parse the GFF annotation file to a hash, also the NEW genbank refseq gff
	# Can extract GFF 'type' such as RNA, gene, CDS, tRNA, rRNA; use regular expression e.g. gene|RNA 	
	# added: aug 2014
	my $gff_filename = shift ;
	my $gff_type = shift ; 
	#print "\tReading gff file $gff_filename\n";
	my @lines = anne_files::read_lines($gff_filename);
	chomp(@lines) ;
	my %result ;
	foreach my $line (@lines) {
		if ($line !~ m/^\#/) {		# line does not start with #
			$line =~ s/\r//g ;
			my @items = split /\t/ , $line ;
			if (@items>2) {
				if ( $items[2] =~ m/($gff_type)/ ) {  # only take this type which match the $gff_type
					my $genome		= $items[0] ;
					my $type		= $items[2] ;
					my $ID	 = '' ;
					if ($items[8] =~ m/ID=(.*?);/)     { $ID = $1 ; }	# gene
					if ($items[8] =~ m/Parent=(.*?);/) { $ID = $1 ; }	# protein
					if ( $ID ne '' ) {
						# 1. parse the columns
						$result{$genome}{$type}{$ID}{genome} 	= $items[0] ;
						$result{$genome}{$type}{$ID}{contig} 	= $items[0] ;
						$result{$genome}{$type}{$ID}{start} 	= int($items[3]) ;
						$result{$genome}{$type}{$ID}{end} 		= int($items[4]) ;
						$result{$genome}{$type}{$ID}{strand} 	= $items[6] ;
						$result{$genome}{$type}{$ID}{length} 	= $result{$genome}{$type}{$ID}{end} - $result{$genome}{$type}{$ID}{start} +1 ;
						$result{$genome}{$type}{$ID}{Name} 			= '' if (!defined($result{$genome}{$type}{$ID}{Name} )) ;
						$result{$genome}{$type}{$ID}{Parent} 		= '' if (!defined($result{$genome}{$type}{$ID}{Parent} )) ;
						$result{$genome}{$type}{$ID}{locus_tag}		= '' if (!defined($result{$genome}{$type}{$ID}{locus_tag} )) ;
						$result{$genome}{$type}{$ID}{old_locus_tag}	= '' if (!defined($result{$genome}{$type}{$ID}{old_locus_tag} )) ;
						$result{$genome}{$type}{$ID}{gene}			= '' if (!defined($result{$genome}{$type}{$ID}{gene} )) ;
						$result{$genome}{$type}{$ID}{product}		= '' if (!defined($result{$genome}{$type}{$ID}{product} )) ;
						$result{$genome}{$type}{$ID}{note}			= '' if (!defined($result{$genome}{$type}{$ID}{note} )) ;
						$result{$genome}{$type}{$ID}{protein_id}	= '' if (!defined($result{$genome}{$type}{$ID}{protein_id} )) ;
						$result{$genome}{$type}{$ID}{geneid}		= '' if (!defined($result{$genome}{$type}{$ID}{geneid} )) ;
						$result{$genome}{$type}{$ID}{gene_biotype}	= '' if (!defined($result{$genome}{$type}{$ID}{gene_biotype} )) ;
						$result{$genome}{$type}{$ID}{pseudo}		= 'false' ;
						# parse the description
						my @descriptions = split ';', $items[8] ;
						foreach my $description (@descriptions) {
							$result{$genome}{$type}{$ID}{Name}			= $1 if ( $description =~ m/^Name=(.*)/ ) ;
							$result{$genome}{$type}{$ID}{Parent} 		= $1 if ( $description =~ m/^Parent=(.*)/ ) ;
							$result{$genome}{$type}{$ID}{locus_tag}		= $1 if ( $description =~ m/^locus_tag=(.*)/ ) ;
							$result{$genome}{$type}{$ID}{old_locus_tag}	= $1 if ( $description =~ m/^old_locus_tag=(.*)/ ) ;
							$result{$genome}{$type}{$ID}{gene} 			= $1 if ( $description =~ m/^gene=(.*)/ ) ;
							$result{$genome}{$type}{$ID}{product}		= $1 if ( $description =~ m/^product=(.*)/ ) ;
							$result{$genome}{$type}{$ID}{note}			= $1 if ( $description =~ m/^Note=(.*)/ ) ;
							$result{$genome}{$type}{$ID}{protein_id}	= $1 if ( $description =~ m/^protein_id=(.*)/ ) ;
							$result{$genome}{$type}{$ID}{gene_biotype}	= $1 if ( $description =~ m/^gene_biotype=(.*)/ ) ;
							$result{$genome}{$type}{$ID}{geneid}		= $1 if ( $description =~ m/^Dbxref=.*GeneID:(.*)/ ) ;
							$result{$genome}{$type}{$ID}{pseudo}		= $1 if ( $description =~ m/^pseudo=(.*)/ ) ;
						}
					}	
				}	
			}	
		}		
	}
	return %result ;
}


sub ptt2hash {
	# store the NCBI ptt file to a hash with the locus tag as key
	# ptt: Location	Strand	Length	PID	Gene	Synonym	Code	COG	Product
	# export: locus(=key), gene, start, end, strand, length, PID, COG, product
	my $pttfile = shift ;
	my %results ;
	my @lines = anne_files::read_lines($pttfile);
	my $tablestart = 0 ;
	foreach my $line (@lines) {
		if ($tablestart) {
			my @items = split /\t/, $line ;
			my @coord = split /\.\./, $items[0] ;
			my $key = $items[5] ; # the locus tag
			$results{$key}{gene} 	= $key ;
			$results{$key}{gene} 	= $items[4] if ($items[4] ne '-');
			$results{$key}{start} 	= $coord[0] ;
			$results{$key}{end} 	= $coord[1] ;
			$results{$key}{strand} 	= $items[1] ;
			$results{$key}{length} 	= $items[2] ;
			$results{$key}{PID} 	= $items[3] ;
			$results{$key}{COG} 	= $items[7] ;
			$results{$key}{product} = $items[8] ;
		}
		if ($line =~ m/^Location/) { $tablestart = 1; }
	}
	return %results ;
}

sub gbk_header2hash {
	my $gbkfile = shift ;
	my %result ;
	my @lines = anne_files::read_lines($gbkfile);
	my @headers = ('LOCUS','DEFINITION','VERSION','ACCESSION','DBLINK','KEYWORDS','SOURCE','ORGANISM');
	foreach my $header (@headers) { $result{$header} = ''; }
	my $getLineage = 0 ;
	my $lineage = '' ;
	foreach my $line (@lines) {
		if ($line =~ m/FEATURES/) { last; }
		foreach my $header (@headers) {				# add the headers to the hash
			if ($line =~ m/.*$header.\s+(.*)/) { 
				$result{$header} = $1 ; 
			}	
		}
		if ($line =~ m/VERSION.*\s+(GI\:.*)/) { # this contains the GI number
			$result{GI} = $1 ;
		}
		if ($line =~ m/REFERENCE/g and $getLineage) {  # the lines between ORGANISM and REFERENCE contains the lineage
			$getLineage = 0 ;
			$lineage =~ s/\s{2,}//g ;
			$result{LINEAGE} = $lineage ;				
			$lineage = '' ;
		}
		if ($getLineage) { $lineage .= $line ; }
		if ($line =~ m/ORGANISM/g)  { $getLineage = 1 ; } # the lines between ORGANISM and REFERENCE contains the lineage
		
	}	
	return %result ;
}

sub get_dna_from_fasta {
	# return a string containing only the DNA
	my $filename = shift ;
	my $result ;
	my @lines = anne_files::read_lines($filename);
	chomp @lines ;
	foreach my $line (@lines) {
		if ($line !~ m/^\>/) { $result .= $line ;	}	
	}
	$result =~ s/(\s+|\r)//g ;
	return $result ;
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
			$key = anne_misc::trim($1) ; # remove leading and taling space
			if ($key =~ /(.*?)\s/) { $key = $1; }
		} else {
			$seq .= $line ;
		}
	}
	# don't forget the last record
	$result{$key} = $seq if ($key ne '') ;
	return %result ;
}

sub g2dfasta2hash {
	# This will parse the fasta file formatted for g2d 
	my $file = shift ;
	my %result ;
	my $seq = '';
	my $key = '' ;
	my $header = '';
	my @lines = anne_files::read_lines($file) ;
	foreach my $line (@lines) {
		if ($line =~ m/^\>(.*)/) { # new record, add the last to the hash
			$result{$key}{seq} = $seq if ($key ne '') ;
			$seq = ''; 
			$header = $1 ;
			my @items = split /\|/ , $header ;
			$key = $items[0] ; 
			$result{$key}{gene} = $items[1] ;
			$result{$key}{proteinid} = $items[2] ;
			$result{$key}{gi} = $items[3] ;
		} else {
			$seq .= $line ;
		}
	}
	# don't forget the last record
	$result{$key}{seq} = $seq if ($key ne '') ;
	my @items = split /\|/ , $header ;
	$key = $items[0] ; 
	$result{$key}{gene} = $items[1] ;
	$result{$key}{proteinid} = $items[2] ;
	$result{$key}{gi} = $items[3] ;
	return %result ;
}

sub get_gbkheader_mysql {
	my $genomename = shift ;
	my %genomeheader ;
	# connect to ProgreSQL for genome data
	my $dbh = DBI->connect("DBI:Pg:dbname=molgenedb;host=localhost", "molgen", "!ggdimme!", {'RaiseError' => 1});
	# Get information on the genome of the input:genomename   e.g. NC_009004
	my $sth = $dbh->prepare("SELECT organism, file, dna_length, cds_count, file FROM dnafeature WHERE accession LIKE '$genomename%'");
	$sth->execute();
	my @vector = $sth->fetchrow ;
	if (!defined($vector[0])) {
		print "ERROR: $genomename not found in the database\n";
		$genomeheader{organism} = '';
	} else {
		print "Organism:$vector[0]\nFile:$vector[1]\nDNA length:$vector[2] bp\n";
		$genomeheader{organism} 	= $vector[0] ;
		$genomeheader{file} 		= $vector[1] ;
		$genomeheader{dna_length} 	= $vector[2] ;
		$genomeheader{cds_count} 	= $vector[3] ;
		$genomeheader{folder} 		= $vector[4] ;
		if ($genomeheader{file} =~ m/(.*)\/.*\.gbk/ ) { $genomeheader{folder} = $1 ; }
	}	
	$sth->finish;
	$dbh->disconnect();
	return %genomeheader ;
}

sub get_genomeheader_from_database {
	my $genomename = shift ;
	my %genomeheader ;
	# connect to ProgreSQL for genome data
	my $dbh = DBI->connect("DBI:Pg:dbname=molgenedb;host=localhost", "molgen", "!ggdimme!", {'RaiseError' => 1});
	# Get information on the genome of the input:genomename   e.g. NC_009004
	my $sth = $dbh->prepare("SELECT organism, file, dna_length, cds_count, file FROM dnafeature WHERE accession LIKE '$genomename%'");
	$sth->execute();
	my @vector = $sth->fetchrow ;
	if (!defined($vector[0])) {
		print "ERROR: $genomename not found in the database\n";
		$genomeheader{organism} = '';
	} else {
		print "Organism:$vector[0]\nFile:$vector[1]\nDNA length:$vector[2] bp\n";
		$genomeheader{organism} 	= $vector[0] ;
		$genomeheader{file} 		= $vector[1] ;
		$genomeheader{dna_length} 	= $vector[2] ;
		$genomeheader{cds_count} 	= $vector[3] ;
		$genomeheader{folder} 		= $vector[4] ;
		if ($genomeheader{file} =~ m/(.*)\/.*\.gbk/ ) { $genomeheader{folder} = $1 ; }
	}	
	$sth->finish;
	$dbh->disconnect();
	return %genomeheader ;
}


sub get_genome_from_database {
	my $genomename = shift ;
	my %genomeheader = get_genomeheader_from_database($genomename);
	my %table ;
	if ($genomeheader{organism} ne '') {

		# connect to ProgreSQL for genome data
		my $dbh = DBI->connect("DBI:Pg:dbname=molgenedb;host=localhost", "molgen", "!ggdimme!", {'RaiseError' => 1});
		# Get information on the genome of the input:genomename   e.g. NC00904
		my $sth = $dbh->prepare("SELECT organism, file, dna_length, cds_count FROM dnafeature WHERE accession LIKE '$genomename%'");
		$sth->execute();
		my @vector = $sth->fetchrow ;



		# Get the annotation of the genes in the genome
		$sth = $dbh->prepare("SELECT locustag, genename, location_start, location_end, strand, protein_id, product, sequence,seq_translation,function,ncbi_geneid FROM genes WHERE organism='".$genomeheader{organism}."'");
		$sth->execute();
		while (my @vector = $sth->fetchrow) {
			$table{$vector[0]}{locus} 	= $vector[0] ;
			$table{$vector[0]}{gene} 	= $vector[1] ;
			$table{$vector[0]}{gene} 	= $vector[0] if ($vector[1] eq "") ;
			$table{$vector[0]}{start} 	= $vector[2] ;
			$table{$vector[0]}{end} 	= $vector[3] ;
			$table{$vector[0]}{strand} 	= $vector[4] ;
			$table{$vector[0]}{protein_id} 	= $vector[5] ;
			$table{$vector[0]}{product} 	= $vector[6] ;
			$table{$vector[0]}{sequence} 	= $vector[7] ;
			$table{$vector[0]}{seq_translation} = $vector[8] ;
			$table{$vector[0]}{function} 		= $vector[9] ;
			$table{$vector[0]}{ncbi_geneid} 	= $vector[10] ;
			$table{$vector[0]}{seq_translation} =~ s/(\{|\})//g ;
			#print "$vector[0]\t$pggenome{$vector[0]}{gene}\t$vector[2]\t$vector[3]\t$vector[4]\n" ;
		}
		$sth->finish;
		$dbh->disconnect();
	}	
	return %table ;
}


sub count_char {
    # returns the number of chars (or regular expression) in the string
	# e.g. for charged residues char = (K|R|H|D|E|C|Y)
	my ($str, $char) = @_;
    my $count = 0 ;
    $count++ while $str =~ /$char/g;
    return $count ;
}

sub calculate_pI {
    # returns the pI of the protein
	my ($seq)  = @_ ;
    my $seqobj = Bio::Seq->new(
				-seq => $seq,
				-id  => 'tocalculate',
				-alphabet => 'protein'
			    ) ; 

    my $calc = Bio::Tools::pICalculator->new(-places => 2);
    $calc->seq( $seqobj );
    my $iep    = $calc->iep ;
    my $charge = $calc->charge_at_pH( 7 );
    return sprintf("%.1f", $iep);;
	#return nearest(.1, $charge), nearest(.1, $iep) ;
}

sub RBS_score {
	# here I expect the small upstream region of the start codon, e.g. 
	my $rbs=shift ;
	my @consensus=split //,"AGGAGG";
	my @seq=split //,$rbs ;
	my $spacing = 5 ; # minimum spacing to the start coding
	my $score = 0 ;
	for my $i (0 .. scalar @seq-6-$spacing) { 
		my $count=0;
		for my $c (0..5) { $count++ if ($seq[$i+$c] eq $consensus[$c]) ; }
		$score=$count if ($count>$score) ;
	}
	return $score ;
} 




sub show_genome {
	foreach my $key (keys %genome) {	
		print "$genome{$key}{gene}\n" if (defined($genome{$key}{gene})) 
    }
}

sub read_condontable {
	open (FILE,'<', $codontablefile) or die ("$codontablefile does not exists") ;
	my @lines = <FILE>;					# Get the codonfile into a array
	chomp(@lines) ; 
	close FILE ;
	foreach my $line (@lines) {
		my @item = split(/\t/, $line) ;
		$codontable{$item[0]} = $item[1];
	}
}

sub translate {
	my $DNA = shift ; 
	$DNA = uc($DNA) ; 		# convert to uppercase
	$DNA =~ s/(\s|\n|\r|\x0d)//g; # remove line breaks spaces etc
	my $protein ;
	my $codon ;
	&read_condontable() ;
	for ( my $i=0; $i<(length($DNA)-2); $i+=3) {
		$codon=substr($DNA,$i,3);
		if (defined $codontable{$codon}) { 
			$protein.=$codontable{$codon}
		}	else {
			$protein.='x' ;
		}
	}
	return $protein ;
}

sub inverse_complement {
	my $DNA = shift ;
	$DNA =~ s/(\s|\n|\r|\x0d)//g;
	my $revcomp = reverse($DNA);
	$revcomp =~ tr/[A,C,G,T,a,c,g,t]/[T,G,C,A,t,g,c,a]/;
	return $revcomp;
}

sub clean_seq {	
	my $DNA = shift ;
	$DNA =~ s/(\s|\x0d|\d)//g;
	return $DNA ;
}	
	
	
sub parsefasta {
	# I do not use BIoPerl over here because I use the preformatted gbk.ffn
	# returns the keys for locus and genenames
	$sessiondir = shift ;
	my $fasta_filename = shift ;
	if (-e $fasta_filename) {
		&add_progress ("Parsing $fasta_filename<br>");
	} else {	
		&add_progress ("File $fasta_filename not found<br>");
		die ;
	}
	# read the file	
	open (FILE,'<', $fasta_filename)  ;
	my @lines = <FILE>;					
	chomp(@lines) ; 
	close FILE ;
	# parse the file
	my $locus ;
	my $genename ;
	foreach my $line (@lines) {
		if ($line =~ m/^\>/) {  # this is the header
			my @header = split(/\|/, $line) ;
			chomp(@header) ; 
			$locus = $header[0] ;
			$locus =~ s/\>// ;  # remove the >
			$genename = $header[1] ;
			$genome{$locus}{gene} 			= $header[1] ;
			$genome{$locus}{proteinID} 		= $header[2] ;
			$genome{$locus}{geneID} 		= $header[3] ;
			$genome{$genename}{locus} 		= $locus ;
			$genome{$genename}{proteinID} 	= $header[2] ;
			$genome{$genename}{geneID} 		= $header[3] ;
		} else {
			$genome{$locus}{sequence} = $line ;
			$genome{$genename}{sequence} = $line ;
		}	
	}
}



sub parsegenbank {
	my $gbk_filename = shift ;
	if (-e $gbk_filename) {
		&add_progress ("Parsing genbank file $gbk_filename<br>");
	} else {	
		&add_progress ("File $gbk_filename not found<br>");
		die ;
	}
    my $bioseq_io = Bio::SeqIO->new( -file=> $gbk_filename) ; 	# read the annotation
	my @ar_pcds  = () ;
	my %prop   = () ;
    my $scaffoldcount = 0 ;
	my $position = 0;
	my $locus ;
	my $taxocount ;

	while( my $seqobj = $bioseq_io->next_seq ) {
		&add_progress("Start scanning the genome<br>");
		$scaffoldcount++ ;
		my @all_fts = $seqobj->get_SeqFeatures ;
		my @tmp_cds = grep { $_->primary_tag eq 'CDS' and $_->has_tag('translation') } @all_fts  ;
		@ar_pcds = () ;
		@ar_pcds   = ( @ar_pcds, @tmp_cds ) ; 
		$genomename{genus}      = $seqobj->species->genus  			if ($seqobj->species->genus) ;
		$genomename{species}    = $seqobj->species->species  		if ($seqobj->species->species) ;
		$genomename{subspecies}	= $seqobj->species->sub_species 	if ($seqobj->species->sub_species) ;
		my @classification = ();
		@classification = $seqobj->species->classification();
		my $taxonomy = '';
		$taxocount = scalar @classification ;
		for (my $i = 1 ; $i <= 6; $i++) {
			$taxonomy .= "$classification[$taxocount-$i]" if (defined($classification[$taxocount-$i])) ;
			$taxonomy .= "\t" ; 
		}
		$genomename{taxonomy} = $taxonomy ;
		if($classification[0] ne '.'){
			$genomename{organism} = $classification[0];
		}else{
			$genomename{organism} = $seqobj->description();
			$genomename{organism} =~ s/|//;
		}
		
		foreach my $pcds ( @ar_pcds ) {
			my %prop = () ;
			$prop{start}  		= $pcds->location->start ; 
			$prop{end}   		= $pcds->location->end ;
			$prop{strand} 		= $pcds->location->strand ;
			$prop{locus}		= ($pcds->get_tag_values('locus_tag'))[0]		if ($pcds->has_tag('locus_tag')) ;
			$prop{gene} 		= ($pcds->get_tag_values('gene'))[0]    		if ($pcds->has_tag('gene')) ;
			$prop{function}		= ($pcds->get_tag_values('function'))[0]     	if ($pcds->has_tag('function')) ;
			$prop{product}		= ($pcds->get_tag_values('product'))[0]     	if ($pcds->has_tag('product')) ;
			$prop{protein_id}	= ($pcds->get_tag_values('protein_id'))[0]  	if ($pcds->has_tag('protein_id')) ;
			$prop{EC_number}	= ($pcds->get_tag_values('EC_number'))[0]  		if ($pcds->has_tag('EC_number')) ;
			$prop{translation}	= ($pcds->get_tag_values('translation'))[0] 	if ($pcds->has_tag('translation')) ;
			$prop{note}			= ($pcds->get_tag_values('note'))[0]        	if ($pcds->has_tag('note')) ;
			$prop{gene_id}	= ( grep{ $_ =~ /^GeneID/ } $pcds->get_tag_values('db_xref') )[0]  if ($pcds->has_tag('db_xref')) ;

			# some GBK don't contain locus info, in that case locus = gene
			$position++ ;
			# in some cases there is no locus or not even a gene defined ==> to solve this:		
			if (!defined $prop{locus}) {
				if (defined $prop{gene}) {
					$prop{locus} = $prop{gene} ;
				} 
				else { # still no locus then give it the protein_id
					if (defined $prop{protein_id}) {
						$prop{locus} = $prop{protein_id} ;
					} 
					else { # unidentified object, just give it a number
						$prop{locus} = 'no_name_'.$position ;
					}	
				} 
			}	
			# handle not defined values
			if (!defined $prop{gene}) 			{ $prop{gene}			= $prop{locus}; }
			if (!defined $prop{function})		{ $prop{function}		='';}  	
			if (!defined $prop{product})		{ $prop{product}		='';}		
			if (!defined $prop{protein_id})		{ $prop{protein_id}	='';}		
			if (!defined $prop{EC_number})		{ $prop{EC_number}='';}	
			if (!defined $prop{translation})	{ $prop{translation}	='';}	
			if (!defined $prop{note})			{ $prop{note}			='';}				
			if (!defined $prop{gene_id})		{ $prop{gene_id}		='';}		
			                                                 
			# add the properties hash to the genome hash
			$locus = $prop{locus} ;
			$genome{$locus} = \%prop ;
			$genome{$locus}{scaffold} = $scaffoldcount ; 
			$genome{$locus}{position} = $position ; 
			$genome{$locus}{context}  = {} ; 	# used later for 3D array to store links to context genes
        }
	}
	&add_progress("Organism= $genomename{organism} => # of scaffolds = $scaffoldcount\n");
}



## the mandatory one (without it no package!!!)
1

