#!/usr/bin/perl

#	 Anne de Jong, October 2016
#		- Get PFAM descriptions
#		

use lib "/usr/molgentools/lib";
use anne_files ;

my $pfam_folder = '/usr/bagel4/Pfam-A';
my $db_name = 'Pfam-A.hmm' ;
my $usage = "Options
	-s	pfam folder [default=$pfam_folder]
	-db	database name [default=$db_name]
	
	e.g. /usr/molgentools/gsea_pro/parse_pfam_descriptions.pl -s /usr/bagel4/Pfam-A -db  Pfam-A.hmm";
	

parseparam() ;

my @lines = anne_files::read_lines("$pfam_folder/$db_name");
my @result = "ACC\tNAME\tDescription";


my $NAME ;
my $ACC ;
foreach my $line (@lines) {
	$NAME = $1 if ($line =~ m/^NAME\s+(.*)/) ;
	$ACC = $1  if ($line =~ m/^ACC\s+(.*)\.\d+/) ;
	if ($line =~ m/^DESC\s+(.*)/) { push @result, $ACC."\t".$NAME."\t".$1;  }
}
print scalar @result ;
print " descriptions found\n";
anne_files::write_lines("$pfam_folder/$db_name.descriptions", @result);


sub parseparam {
    my $var ;
    my @arg = @ARGV ;
    while(@arg) {
        $var = shift(@arg) ;
        die $usage if ($var eq '-h' or $var eq '--help') ;
		$pfam_folder	= shift(@arg) if($var eq '-s') ;
		$db_name		= shift(@arg) if($var eq '-db') ;
    }
}

