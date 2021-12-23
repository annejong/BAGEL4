#!/usr/bin/env perl

use strict ;
use warnings ;
use lib "/data/bagel4/lib" ;
use bagel4_functions ;
use lib "/data/molgentools/lib";
use anne_files ;
use anne_genomics ;
use File::Basename;
use LWP::Simple;


#my $xml = get_uniprot_xml('LANN_LACLA') ;
#anne_files::write_string('LANN_LACLA.xml', $xml) ;

my $xml = join "\n", anne_files::read_lines('LANN_LACLA.xml') ;
#print anne_files::print_lines(@xml) ;
#print $xml ;

#while ($xml =~ m/\<feature type=(.*)\<\/feature\>/g) {
$xml =~ s/\n//g ;
my $index = 0 ;
while ((substr $xml, $index) =~ m/<feature (.*?)<\/feature>/) {
	my $feature = $1 ;
	$index += $-[0] + 1 ;
	my $type = ''; 			$type = $1 			if ($feature =~ m/type=\"(.*?)\"/) ;
	my $description = ''; 	$description = $1	if ($feature =~ m/description=\"(.*?)\"/) ;
	my $beginpos = ''; 	$beginpos = $1.' -'	if ($feature =~ m/begin position=\"(.*?)\"/) ;
	my $endpos = ''; 	$endpos = $1	if ($feature =~ m/end position=\"(.*?)\"/) ;
	my $position = ''; 	$position = $1	if ($feature =~ m/position position=\"(.*?)\"/) ;
	print "$type $description $position $beginpos $endpos\n";
}	

sub get_uniprot_xml {
	my $uniprot = shift ;
	my $url = "http://www.uniprot.org/uniprot/$uniprot.xml" ;
	my $content = get($url);
	return $content ;
}
