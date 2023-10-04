package anne_files ;

##############################################################

## 	Anne de Jong
##	University of Groningen
##	the Netherlands
##	anne.de.jong@rug.nl

##############################################################
##
##   
##  2013 januari,  module for file handling
##

use strict ;
use warnings ;
use File::Find qw(find);



BEGIN {
    use Exporter ();
    @kleur::ISA         = qw(Exporter);
    @kleur::EXPORT      = qw();
}


sub in_array {
	# check if $str is in array
	my ($str, @array) = @_ ;
	my $result = 0 ;
	$result = 1 if ( grep( /^$str$/, @array ) ) ;
	return $result ;
}


sub dec2bin {
   my $str = unpack("B32", pack("N", shift));
   $str =~ s/^0+(?=\d)//;   # otherwise you'll get leading zeros
   return $str;
}

sub unique_array {
	# remove duplicates / replicates from array
    return keys %{{ map { $_ => 1 } @_ }};
}

sub unique_pipe_list {
	#split the string, make unique, and join ;
	return join "|", unique_array(split /\|/, shift)  ;
}

sub chk_php_code {
	# security routine to check for PHP code
	my $filename = shift ;
	my $result = 0 ;
	open(FILE, $filename) or die "Cannot find file\n" ;
	my $content = <FILE>;
	if ($content =~ /\<\?php/) {
		$result = 1 ;
		last ;
	}	
	close(FILE);
	return $result ;
}	

sub reformat_json {
	# make it more human readable
	my $inputfilename = shift ;	
	my @lines = anne_files::read_lines($inputfilename) ;
	$lines[0] =~ s/\{/\{\n/g ;
	$lines[0] =~ s/\}/\n\}/g ;
	$lines[0] =~ s/\,/\,\n/g ;
	anne_files::write_lines($inputfilename,@lines) ;
}


sub table2html {
	my $inputfilename = shift ;
	my $outputfilename = shift ;
	if (!-e $inputfilename) { 
		print "WARNING: $inputfilename not found\n";  
	} else {
		my @lines = read_lines($inputfilename) ;
		my @html ;
		my $color_odd	= "#B7CEEB" ;
		my $color_even	= "#ADD8E6" ;
		my $odd = 0 ;
		my $count = 1 ;
		push @html, "<TABLE>\n" ;
		foreach my $line (@lines) {
			my @elements = split /\t/, $line ;
			if ($count == 1) { # print header
				push @html, "<thead>\n\t<tr bgcolor=#4682B4>" ;
				foreach (@elements) { push @html, "<th><font color=#E0FFFF>$_</font></th>"; }
				push @html, "</tr>\n</thead><tbody>\n" ;
			} else {
				my $bg = $color_odd ; 
				if ($odd) { $odd= 0; $bg = $color_even; } else { $odd= 1; }
				push @html, "\t<tr bgcolor=$bg>" ;
				foreach (@elements) { push @html, "<td style=padding-left:2px;padding-right:15px;>$_</td>"; }
				push @html, "</tr>\n" ;
			}	
			$count++ ;
		}	
		push @html, "</tbody></TABLE>\n" ;
		write_lines($outputfilename, @html) ;
	}
}



sub get_files_from_dir_regex {
	# regex added to get_files_from_dir
	my ($directory, $regex) = @_ ;
	opendir(DIR,$directory);
	my @files = readdir(DIR);
	closedir(DIR);
	my @results ;
	foreach(@files){ push @results, $_ if ($_ =~ m/$regex/) ; }	
	return @results ;
}

sub get_files_from_dir {
	# OLD version without regex
	my $directory = shift ;
	my @results ;
	opendir (DIR, $directory) or die $!;
    while (my $file = readdir(DIR)) {
        # We only want files
        next unless (-f "$directory/$file");
        # Use a regular expression to ignore files beginning with a period
        next if ($file =~ m/^\./);
		push @results, $file ;
    }
	closedir(DIR);
	return @results ;
}

sub get_files_from_subdirs {	# searches the path and its sub directories for files matching with $regexp
	my ($path, $regexp) = @_ ;									# parse the commandline parameters	
	print "====> $path\n";
	print "====> $regexp\n";
	my $files = [] ;
	($regexp, $path, $files) = &filesearch($regexp, $path, $files) ;
	return( $files ) ;											# return the previous info
}


sub get_files_from_subdirs_v2 {	# NEW July 2017: return files in subdirs
	my $folder = shift;
	my $regex = shift ;
	my @files;
	push @files, list_dir($folder, $regex);
	return @files ;
}

sub list_dir {
    #my @dirs = @_;
    my $dir = shift ;
	my $regex = shift ;
	my @files;
    find({ wanted => sub { push @files, glob "\"$_\"" } , no_chdir => 1 }, $dir);
	my @result ;
	#foreach my $file (@files) { push @result, $file if ($file =~ m/$regex/) ; }
	#foreach my $key (@files) { print "fasta file=$key\n"; }
    return @files;
}

sub filesearch {														# use an anonymous function to limit code complexity
	my ($regexp, $path, $files) = @_ ;									# parse the parameters							
	$path .= '/' if ($path !~ /\/$/) ;									# add a slash to the end of the path if it is not provided
	opendir( my $ha, $path ) ;											# open the directory
	my @all = sort map {$path.$_} grep {$_ !~ /^\./} readdir($ha) ;	# get all the files and subdirectories in the directory 
	closedir( $ha )	;													# close the directory
	my @tocheck = grep {-d $_ } @all ;									# get the subdirectories which need checking	
	my @toadd   = grep {$_ =~ /$regexp/} @all ;							# get the files to add 	
	foreach my $fn (@toadd ) {
		push @$files, $fn ;												# add the files to the store
	}		
	foreach my $sdir ( @tocheck ) {										# check the subdirectories	
		($regexp, $path, $files) = &filesearch($regexp, $sdir, $files) ;
	}
	return($regexp, $path, $files) ;									# return the results
} 

sub get_filename_from_path {
	# return the filename from a full path string
	my $file = shift ;
	my @items = split (/\//, $file) ;
	my $result = pop(@items) ;
	return $result ;
}

sub get_folder_from_path {
	# return the folder from a full path string
	my $file = shift ;
	my @items = split (/\//, $file) ;
	pop(@items) ;
	my $result = join ("/", @items) ;
	return $result ;
}	


sub write_string {
	my ($filename, $string) = @_;
	open (FILE,">$filename") or die ("Could not write to $filename");
	print FILE $string;
	close FILE ;
}

sub write_lines {
	my ($filename, @lines) = @_ ;
	chomp @lines ;
	open (FILE,">$filename") or die ("Could not write to $filename");
	foreach my $line (@lines) { 
		$line =~ s/\r//g ;
		print FILE $line."\n" ; 
	}
	close FILE ;
}

sub append_lines {
	my ($filename, @lines) = @_ ;
	open (FILE,">>$filename") or die ("Could not write to $filename");
	foreach my $line (@lines) { print FILE $line."\n" ; }
	close FILE ;
}

sub read_lines {
	my $filename=shift;
	my @lines = ''; 
	if (-e $filename) {
		open (FILE,"<$filename") or die ("Could not read $filename");
		@lines = <FILE>;	
		close FILE ;
		chomp @lines ;
	}
	#for( my $i = 0; $i < scalar @lines; $i++ ){ $lines[$i] =~ s/(\n|\r|\x0d)//g; }		 # remove line breaks etc
	return @lines ;
}

sub show_lines {
	my @lines = @_ ;
	foreach my $line (@lines) { print $line."\n"; }
}


sub properfilename {
	my $str = shift ;
	# these chars are not allowed in filenames - \ / | : , ; * " ? < > ( )
	$str =~ s/(\-|\ |\,|\s|\t|\\|\/|\||\:|\;|\*|\"|\?|\<|\>)/\_/g ; 
	$str =~ s/(\(|\))//g ;
	$str =~ s/\.//g ;
	$str =~ s/\.$//g ;
	$str = substr ($str, 1, 200) if (length($str) > 200) ; 
	return $str;
}	

sub dos2linux {
	# remove ^M etc.. 
	my $filename = shift ;
	my @lines = read_lines($filename);
	chomp @lines ;
	write_lines($filename,@lines) ;
}

sub mac2linux {
	my $filename = shift ;
	my @lines = read_lines($filename);
	my @result ;		
	foreach my $line (@lines) {
		$line =~ tr/\cM/\n/ ;
		push @result, $line;
	}	
	write_lines($filename, @result) ;
}




sub cleanFile {
	# remove ^M etc.. 
	my $filename = shift ;
	my @lines = read_lines($filename);
	chomp @lines ; # remove ^M chars 
	foreach my $line (@lines) {
		$line =~ s/^\s+|\s+$//g ; #  trim space at ends
	}
	write_lines($filename,@lines) ;
}


sub add_header {
	# add the string to the top of the file, only if the header does not exists
	my $filename = shift ;
	my $str = shift ;
	my @lines =  read_lines($filename) ;
	@lines = $str if (!defined($lines[0]));
	if ($str ne $lines[0]) {
		unshift (@lines, $str."\n" ) ;
		write_lines ($filename, @lines) ;
	}	
}

sub add_keys {
	# If a file does not contain unique keys, this routine generate a key column
	my $filename = shift ;
	my @lines =  read_lines($filename) ;
	my $count = 0 ;
	my @newlines ;
	foreach my $line (@lines) {
		push @newlines, "KEY_$count\t$line\n" ;
		$count++;
	}
	write_lines($filename, @newlines) ;
}

sub get_table_header {
	my $filename = shift ;
	open( my $file, $filename);
	my @header = split '\t', <$file>;
	chomp @header ;
	#my @lines =  read_lines($filename) ;
	#chomp @lines ;
	#my @tmp = split "\t", $lines[0];
	return @header ;
}



sub Table2hash_v2 { 
	# new version. Allows adding rowID
	# generic routine to fill a hash from a table file without ID
	# the file should be tab delimited and contain a header. The first columns contains the keys of the hash if $addID=false
	my $filename = shift ;
	my $addID = shift ;  # true|false
	my @lines = anne_files::read_lines($filename) ;
	my %result ;
	$lines[0] =~ s/\ //g ; # remove spaces from the header
	$lines[0] =~ s/\r//g ; # remove return
	my @header = split /\t/, shift @lines ;
	my $rowcount = 0 ;
	foreach my $line (@lines) {
		my @col = split /\t/, $line ;
		if (scalar @col > 1 and $col[0] ne '') { 
			$rowcount++ ;
			my $col_count=0;
			foreach my $columname (@header) {
				my $key = $col[0] ;
				$key = $rowcount if ($addID eq 'true') ;
				$result{$key}{$columname} = $col[$col_count]; 
				$col_count++;
			}
		}	
	}
	return %result ;
}



sub read_table_to_hash { 
# generic routine to fill a hash from a table file 
# the file should be tab delimited and contain a header. The first columns contains the keys of the hash
	my $filename = shift ;
	my %tmp ;
	open(FILE, "$filename" ) or die("could not find the file $filename\n") ;
	my(@lines) = <FILE>;
	chomp @lines ;
	@lines = '' if (!defined($lines[0]));
	$lines[0] =~ s/\ //g ; # remove spaces from the header
	$lines[0] =~ s/\r//g ; # remove return
	my @header = split /\t/, $lines[0] ;
	#foreach my $columname (@header) { print "$columname\n"; }
	my $headerline = 1 ;
	foreach my $line (@lines) {
		if ($headerline) {
			$headerline = 0 
		} else {
			$line =~ s/\r//g ;
			my @col = split /\t/, $line ;
			if (scalar @col > 1 and $col[0] ne '') { 
				my $col_count=0;
				foreach my $columname (@header) {
					$tmp{$col[0]}{$columname} = $col[$col_count]; 
					$col_count++;
				}
			}	
		}
	}
	close(FILE) ;
	return %tmp ;
}

sub read_plain_table_to_hash {
	# no header and no key, just the columns as numbers
	my $filename = shift ;
	my @lines = read_lines($filename) ;
	my $key = 0 ;
	my %results ;
	foreach my $line (@lines) {
		my @items = split "\t", $line ;
		for (my $col = 0 ; $col<scalar(@items); $col++ ) { $results{$key}{$col} = $items[$col] ; 	}	
		$key++ ;
	}	
	return %results ;	
}

sub get_table_header_from_hash {
	my %table = @_ ;
	my $first = (keys %table)[0] ;
	my @result = (keys %{$table{$first}}) ;
	return @result ;
}

sub print_table_of_hash_of_hashes {
	my %table = @_ ;
	my @headers = sort (get_table_header_from_hash(%table));
	print "key\t".join("\t",@headers)."\n";
	foreach my $key (sort keys %table )  {
		my @line ;
		foreach my $header (@headers) {
			push @line, $table{$key}{$header};
		}
		print "$key\t".join("\t",@line)."\n";
	}
}

sub print_lines {
	my @lines = @_ ;
	foreach my $line (@lines) { print $line."\n"; }
}

sub save_table_of_hash_of_hashes {
	# returns the lines to be saved
	my %table = @_ ;
	my @result ;
	my @headers = sort (get_table_header_from_hash(%table));
	push @result, "key\t".join("\t",@headers);
	foreach my $key (sort keys %table )  {
		my @line ;
		foreach my $header (@headers) {
			push @line, $table{$key}{$header};
		}
		push @result, "$key\t".join("\t",@line);
	}
	return @result ;
}

sub print_hash_of_hashes {
	my %table = @_ ;
	foreach my $key (sort keys %table )  {
		my @line ;
		foreach my $key2 (sort keys %{$table{$key}}) {
			push @line, "$key2=$table{$key}{$key2}";
		}
		print "$key\t".join("\t",@line)."\n"; ;
	}
}


sub get_second_hashkey {
	my %table = @_ ;
	my @keys = keys %table ;
	my @keys2 = (keys %{$table{$keys[1]}}) ;
	return @keys2 ;
}


sub write_log {
	my $filename=shift;
	my $string=shift;
	my $print=shift;	# true means print to screen
	open (FILE,">>$filename") or die ("Could not write to $filename");
	print FILE "$string\n";
	close(FILE);
	print "$string\n" if ($print eq 'true');
	my $tmp = 1;
	return $tmp ;
}



sub old_write_lines_to_htmltable {
	my ($filename, @lines) = @_ ;
	chomp @lines ;
	my @elements ;
	my $count = 1 ;
	my @html = "	<style>
	#g2dTable1 {
		padding-bottom:1em;
		border-spacing:.5rem;
	}
	th {
		background-color:#2e4f17;
		text-align:left;
		color:#d5e1cc;
	}
	td {
		background-color:#bad8a4;
		text-align:right;
		color:#002927;
		padding:.5rem;
	}
	</style>
	<TABLE id=g2dTable1>\n";
	#my @html = "<table>\n" ;
	foreach my $line (@lines) {
		@elements = split /\t/, $line ;
		if ($count == 1) { # print header
			push @html, "<thead>\n\t<tr>" ;
			foreach (@elements) { push @html,  "<th>$_</th>"; }
			push @html,  "</tr>\n</thead><tbody>\n" ;
		} else {
			push @html,  "\t<tr>" ;
			foreach (@elements) { push @html,  "<td>$_</td>"; }
			push @html,  "</tr>\n" ;
		}	
		$count++ ;
	}	
	push @html,  "</tbody></table>\n" ;
	write_lines($filename, @html) ;
}

sub write_lines_to_htmltable {
	my ($filename, @lines) = @_ ;
	chomp @lines ;
	my $style = "<style>
		#g2dTable1 {
			padding-bottom:1em;
			border-spacing:.5rem; }
		th {background-color:#2e4f17;
			text-align:left;
			color:#d5e1cc; }
		td {background-color:#bad8a4;
			text-align:right;
			color:#002927;
			padding:.5rem; }
	</style>";
	my @html = "$style\n <TABLE id=g2dTable1>\n";
	# the header
	my @elements = split /\t/, shift @lines ;
	push @html, "<thead>\n\t<tr>" ;
	foreach (@elements) { push @html,  "<th>$_</th>"; }
	push @html,  "</tr>\n</thead><tbody>\n" ;
	# the body
	foreach my $line (sort @lines) {
		@elements = split /\t/, $line ;
		push @html,  "\t<tr>" ;
		foreach (@elements) { push @html,  "<td>$_</td>"; }
		push @html,  "</tr>\n" ;
	}	
	push @html,  "</tbody></table>\n" ;
	write_lines($filename, @html) ;
}



## the mandatory one (without it no package!!!)
1