#!/usr/bin/env perl

# Anne de Jong, October 2017
# Add promoters and terminiators to json

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
my $queryname ;


my $usage = "option:
		-s		sessiondir [default=$sessiondir]
		-queryname queryname of the table names. e.g. AOI_1
			
		reads the promoters and terminators table and convert it to a JSON file	
			promoters tab delimited table [pos strand description]  queryname.promoters
			terminators tab delimited table [pos strand description] queryname.transterm.table
		
		e.g.  /usr/bagel4/bagel4_promoter_term_2_json.pl -s /usr/bagel4/test -queryname NC_0085331
" ;

my $yline_hight = 10 ;
my $prom_hight = 2 ;
my $prom_point = 0.4 ;

my $term_hight = 2 ;
my $term_loop = 0.4 ;

&parseparam() ;

# ----------------------------------------------------------- main ------------------------------------------------------------------------------

print "Adding Promoter and Terminator data to JSON\n";

	
my %fasta = anne_genomics::fasta2hash("$sessiondir/$queryname.AOI.fna") ;  # read the MULTI fasta file AND the AOI names

my $yline = 10;
foreach my $queryname (sort keys %fasta) {
	print "\tAdding $queryname\n";
	my $region_size = length($fasta{$queryname}) ;  # get the region size for scaling
	my @promoter_elements ;
	my @terminator_elements ;
	# Add the promoters elements
		my @lines = anne_files::read_lines("$sessiondir/$queryname.promoters");
		foreach my $line (@lines) {
			my @items = split "\t", $line ;
			if (scalar @items >2) {	  
				print "\t\tpromoter\t$line\n";
				my $position = 100 * $items[0] / $region_size  ;  # rescale to 100
				my @element = '{' ;
				push @element, '"name": "promoter",' ;
				push @element, '"color": "#2d65bf",' ;
				push @element, '"yline": '.$yline.',';
				push @element, "\"xtext\": $position,";
				push @element, "\"annotation\":\"$items[2]\",";
				push @element, "\"points\": [".add_promoter_polygon($position, $items[1])."]}" ;
				push @promoter_elements, (join "\n\t", @element)."\n" ;  # add each element to elements
			}	
		}


	# Add the terminators elements
		@lines = anne_files::read_lines("$sessiondir/$queryname.transterm");
		foreach my $line (@lines) {
			if ($line =~ m/^(\d+)\t(.)\t(.*)/) {
				print "\t\tterminator\t$line\n";
				my $position = 100 * $1 / $region_size  ;  # rescale to 100
				my $ori = 1; $ori = -1 if ($2 eq '+') ;
				my @element = '{' ;
				push @element, '"name": "terminator",' ;
				push @element, '"color": "#ba4a28",' ;
				push @element, '"yline": '.$yline.',';
				push @element, "\"xtext\": $position,";
				push @element, "\"annotation\":\"$3\",";
				push @element, "\"stem\": [".json_coord($position,0).",".json_coord($position,$term_hight*$ori)."]," ;
				push @element, "\"loop\": ".json_circle($position,$term_hight*$ori,$term_loop)."}" ;
				push @terminator_elements, (join "\n\t", @element)."\n" ;  # add each element to elements
			}	
		}
	# make the final JSON	
		my @json = '{ '; # json start and end with { }
		push @json, '"promoters": [';
			push @json, join ',', @promoter_elements ;
			push @json, "]," ;
		push @json, '"terminators": [';
			push @json, join ',', @terminator_elements ;
			push @json, "]";
		push @json, "\n}" ;  # end JSON with }

	anne_files::write_lines("$sessiondir/$queryname.prom_term.json", @json) ;
}	

#print (join /\t/,@json) ;
#print "\n";
# ----------------------------------------------------------------------------- functions ----------------------------------------------------------


sub json_circle {
	my ($cx,$cy, $r) = @_ ;
	return "{\"cx\":$cx,\"cy\":$cy,\"r\":$r}";
}


sub json_coord {
	my ($x,$y) = @_ ;
	return "{\"x\":$x,\"y\":$y}";
}

sub add_promoter_polygon {
	my ($position, $strand)  = @_ ;
	my @result ;
	if ($strand eq '+') {
		push @result, json_coord($position,0) ;
		push @result, json_coord($position,-$prom_hight) ;
		push @result, json_coord($position+$prom_point, -$prom_hight+($prom_point/2)) ;
		push @result, json_coord($position, -$prom_hight+$prom_point) ;
	} else {
		push @result, json_coord($position,0) ;
		push @result, json_coord($position,$prom_hight) ;
		push @result, json_coord($position-$prom_point, $prom_hight-($prom_point/2)) ;
		push @result, json_coord($position, $prom_hight-$prom_point) ;
	}
	return join ',', @result ;
}



sub parseparam {
    my $var ;
    my @arg = @ARGV ;
    while(@arg) {
        $var = shift(@arg) ;
        die $usage if ($var eq '-h' or $var eq '--help') ;
		$sessiondir		= shift(@arg) if($var eq '-s') ;
		$queryname		= shift(@arg) if($var eq '-queryname') ;
    }
	die $usage if (!$queryname) ;
}


