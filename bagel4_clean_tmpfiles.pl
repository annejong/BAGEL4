#!/usr/bin/env perl

use strict ;
use warnings ;
use lib "/data/molgentools/lib";
use anne_files ;


my $sessiondir ;

my $usage = "options:
				-s Sessiondir 
		e.g.  /data/bagel4/bagel4_clean_tmpfiles.pl -s /tmpdrive/bagel4wrapper/129.125.143.3510il8du3k1lhd14gq06to6odf0500 
		
		remove empty files
		
\n - Anne de Jong - \n" ;



&parseparam() ;


my @filenames = anne_files::get_files_from_subdirs_v2($sessiondir, '*') ;
my $total = scalar @filenames ;
my $removed = 0 ;
foreach my $file (@filenames) {
	my $size = -s $file ;
	my $del = 0 ;
	$del = 1 if ($size < 1350 and $file =~ m/domtblout$/) ;
	$del = 1 if ($size < 50 and $file =~ m/html$/) ;
	$del = 1 if ($file =~ m/\.(phr|pin|psd|psi|psq|six_frames)/) ;
	$del = 1 if ($size < 1) ;
	if ($del) {
		print "$file\t$size\n";
		unlink $file ;
		$removed++;
	}	
}

print "\n\nTotal files= $total, removed= $removed\n";

sub parseparam {
    my $var ;
    my @arg = @ARGV ;
    while(@arg) {
        $var = shift(@arg) ;
        die $usage if ($var eq '-h' or $var eq '--help') ;
        $sessiondir 	= shift(@arg) if($var eq '-s') ;
    }
    die $usage if (!$sessiondir) ;
}

