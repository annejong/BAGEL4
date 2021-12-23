#!/usr/bin/env perl


# Convert the bagel4.conf to json

# depends on MOODS 


use strict ;
use warnings ;
use lib "/data/bagel4/lib" ;
use bagel4_functions ;
use lib "/data/molgentools/lib";
use anne_files ;
use anne_misc ;
use File::Basename;


# ---------------------------------------------------------- parameters -------------------------------------------------------------------------

my $sessiondir = '.';
my $program_dir = dirname($0) ;
my $fna_file  ;
my $outfile = 'promoterscan.table' ;
my $usage = "option:
		-s		sessiondir [default=$sessiondir]

				
	Conver the bagel4.conf to bagel4.conf.json
				
	e.g.  /data/bagel4/bagel4_conf_2_json.pl -s /tmp/BAGEL4WRAPPER/test	
" ;
my %conf = bagel4_functions::read_conf("$program_dir/bagel4.conf") ;

&parseparam() ;

# ---------------------------------------------------------- main -------------------------------------------------------------------------


my @fields ;
foreach my $key (keys %conf) {
	my $c = '"';
	if (anne_misc::is_numeric($conf{$key})) { $c = ''; }
	push @fields, "\t\"$key\": $c$conf{$key}$c";
}	
my @json = '{"bagel4_conf": {';
push @json, join ",\n", @fields ;
push @json, "}}";
	
anne_files::write_lines("$sessiondir/bagel4.conf.json",@json); 	

# ---------------------------------------------------------- functions -------------------------------------------------------------------------

sub parseparam {
    my $var ;
    my @arg = @ARGV ;
    while(@arg) {
        $var = shift(@arg) ;
        die $usage if ($var eq '-h' or $var eq '--help') ;
		$sessiondir		= shift(@arg) if($var eq '-s') ;
    }
}

