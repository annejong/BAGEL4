#!/usr/bin/env perl

#  BAGEL4 program 
# 	Anne de Jong and Auke van Heel
#	University of Groningen
#	the Netherlands
#	anne.de.jong@rug.nl
#
#  	BedGraph 2 json
#
# v 1.0 
#
#  15 Feb 2018: Some bed graph files do not contain a start but only end, so we use the 3th column 


use strict ;
use warnings ;
use lib "/data/molgentools/lib";
use anne_files ;
use anne_genomics ;



# -----------------------------------------------------------------------------------  global vars ---------------------------------------------------------------------------
my $sessiondir = "." ;
my $bedgraph ;
my $queryname = 'my_name';
my $start = 100 ;
my $end = 500 ;
my $outputfile = 'my_bedgraph.json' ;
my $usage = "bagel4_bedgraph_2_json.pl
				-s Sessiondir [default=$sessiondir]
				-bed full filename of BedGraph file
				-queryname queryname
				-start start position 
				-end end position
				-o outputfilename
				
e.g.  /data/bagel4/bagel4_bedgraph_2_json.pl -bed /data/www/jbrowse/Tonia_Sporulation/data/dataset02/4140_T2A.plus.sort.bedgraph -start 60 -end 240 
	
" ;

&parseparam() ;


# ---------------------------------------------------------- main -------------------------------------------------------------------------

# 1. Read the BedGraph
	my @lines = anne_files::read_lines($bedgraph);

# 2. Get range of the BedGraph
	my $range = $end - $start ;
	my %bedgraph ;  # start positions are the keys
	foreach my $line (@lines) {
		my @items = split "\t", $line ;
		my $position = $items[2] ;
		if (scalar @items >2) {	
			if ($position > $start and $position < $end) { 
				$bedgraph{$position} = abs($items[3]) ; 
			}
		}	
	}
	my @values = sort { $bedgraph{$a} <=> $bedgraph{$b} } keys %bedgraph;
	my $maxy_log = log(1+$values[-1]) ;  # the last key
	print "range = $range; maxSignal = $values[-1]\n" ;
		

# 3. Scale the points	
	my @points = '{"x":'.ScaleX($start).',"y":'.ScaleY(0).'}' ;  # the start point
	foreach my $pos (sort keys %bedgraph) {
		push @points, '{"x":'.ScaleX($pos).',"y":'.ScaleY($bedgraph{$pos}).'}' ;
	}
	push @points, '{"x":'.ScaleX($end).',"y":'.ScaleY(0).'}' ;


# 4. Add the points to JSON
	my $all_points = join ",",@points ;
	my @json = '{' ;
	push @json, '"RNA": [{ ' ;
	push @json, '"name": "mRNA",' ;
	push @json, '"color": "#2d65bf",' ;
	push @json, '"points":['.$all_points.']}]' ;
	push @json, '}' ;

# 5. write the json
	anne_files::write_lines("$sessiondir/$queryname.bedgraph.json", @json);

# ---------------------------------------------------------- functions -------------------------------------------------------------------------
	
sub ScaleX {
	my $val = shift ;
	return 100 * ($val-$start)/$range ;
}	

sub ScaleY {
	my $val = abs(shift) ;
	$val = log(1+$val) ;
	return 20 * $val/$maxy_log ;
}	
	
sub parseparam {
    my $var ;
    my @arg = @ARGV ;
    while(@arg) {
        $var = shift(@arg) ;
        die $usage if ($var eq '-h' or $var eq '--help') ;
		$sessiondir		= shift(@arg) if($var eq '-s') ;
		$bedgraph 		= shift(@arg) if($var eq '-bed') ;
		$queryname		= shift(@arg) if($var eq '-queryname') ;
		$start		= shift(@arg) if($var eq '-start') ;
		$end		= shift(@arg) if($var eq '-end') ;
    }
	die $usage if (!$bedgraph) ;
}
