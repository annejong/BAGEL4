package bagel4_blast ;

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
##  2017 July
##  Package for parsing blast results for BAGEL4
##

use strict ;
use warnings ;
use lib "/usr/molgentools/lib";
use anne_files;
use LWP::Simple;

BEGIN {
    use Exporter ();
}



# ------------------------------------------------------------------  functions  ----------------------------------------------------------------


sub get_uniprot_features {
	my $uniprotXMLfile = shift ;  # full path to the uniprot xml filename
	
	my @lines = anne_files::read_lines($uniprotXMLfile) ;
	my $xml = join '',@lines ;
	$xml =~ s/\n//g ; # remove all returns 
	my $index = 0 ;
	my @results ;
	my %features ;
	my $key = 0 ;
	while ((substr $xml, $index) =~ m/<feature (.*?)<\/feature>/) {
		my $feature = $1 ;
		$index += $-[0] + 1 ;
		$key++ ;
		my $type = ''; 			$type = $1 			if ($feature =~ m/type=\"(.*?)\"/) ;
		my $description = ''; 	$description = $1	if ($feature =~ m/description=\"(.*?)\"/) ;
		my $beginpos = ''; 	$beginpos = $1	if ($feature =~ m/begin position=\"(.*?)\"/) ;
		my $endpos = ''; 	$endpos = $1	if ($feature =~ m/end position=\"(.*?)\"/) ;
		my $position = ''; 	$position = $1	if ($feature =~ m/position position=\"(.*?)\"/) ;
		push @results, "$type $description $position $beginpos $endpos";
		$features{$key}{type} = $type ;
		$features{$key}{description} = $description ;
		$features{$key}{position} = $position ;
		$features{$key}{beginpos} = $beginpos ;
		$features{$key}{endpos} = $endpos ;
	}		
	return %features ;
}



sub get_uniprot {
	my $uniprot = shift ;
	my $url = "http://www.uniprot.org/uniprot/$uniprot.xml" ;
	my $xml = get($url);
	$xml =~ s/\n//g ; # remove all returns 
	my $index = 0 ;
	my @results ;
	my %features ;
	my $key = 0 ;
	while ((substr $xml, $index) =~ m/<feature (.*?)<\/feature>/) {
		my $feature = $1 ;
		$index += $-[0] + 1 ;
		$key++ ;
		my $type = ''; 			$type = $1 			if ($feature =~ m/type=\"(.*?)\"/) ;
		my $description = ''; 	$description = $1	if ($feature =~ m/description=\"(.*?)\"/) ;
		my $beginpos = ''; 	$beginpos = $1.' -'	if ($feature =~ m/begin position=\"(.*?)\"/) ;
		my $endpos = ''; 	$endpos = $1	if ($feature =~ m/end position=\"(.*?)\"/) ;
		my $position = ''; 	$position = $1	if ($feature =~ m/position position=\"(.*?)\"/) ;
		push @results, "$type $description $position $beginpos $endpos";
		$features{$key}{type} = $type ;
		$features{$key}{description} = $description ;
		$features{$key}{position} = $position ;
		$features{$key}{beginpos} = $beginpos ;
		$features{$key}{endpos} = $endpos ;
	}		
	return %features ;
}




## the mandatory one (without it no package!!!)
1
