package bagel4_functions ;

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
##  Package with new functions for BAGEL4
##

use strict ;
use warnings ;
use lib "/usr/molgentools/lib";
use anne_files;
use anne_genomics;
use anne_misc;
use LWP::Simple;

BEGIN {
    use Exporter ();
}

my %codontable ;
my $codontablefile = '/usr/bagel4/tables/codon_table.txt';
my $logfile 		= "logfile_bagel.txt";
my $sessionprogress = "sessionprogress.txt";


# ------------------------------------------------------------------  functions  ----------------------------------------------------------------


sub get_sessionID {
	# return only the sessionID from the sessiondir
	my $sessiondir = shift ;
	my @items = split /\//,$sessiondir ;
	return $items[-1] ;
}


sub write_bagel_log {
	# write to log file, sceen, and web-sessionprogress
	# e.g. anne_files::write_bagel_log( "ORF calling on AOI $key", $sessiondir,,1,1,1 );
	my ($my_str, $sessiondir, $screen, $log, $sessionprogress) = @_ ;
	print "$my_str\n" if ($screen) ;
	anne_files::append_lines("$sessiondir/$logfile", $my_str ) if ($log) ;
	anne_files::append_lines("$sessiondir/sessionprogress"), $my_str if ($sessionprogress) ;
}

sub read_conf {
	# returns the configuration file as a hash table
	my $filename = shift ;
	my %results ;
	open(FILE,"$filename") or die('could not find configuration file $filename\n') ;
	my(@lines) = <FILE>;
	chomp @lines ;
	foreach my $line (@lines) {
		$line =~ s/(\s|\t)//g ; # remove spaces and tabs
		if ($line =~ m/(.*)\=(.*)/g ) {	$results{$1}=$2; }	
	}
	return %results ;
}


sub read_condontable {
	my @lines = anne_files::read_lines($codontablefile) ;
	foreach my $line (@lines) {
		my @item = split(/\t/, $line) ;
		$codontable{$item[0]} = $item[1];
	}
}


sub six_frames {
	# translate DNA to 6 frames
	my $DNA = uc(shift) ;
	$DNA =~ s/(\s|\n|\r|\x0d)//g; # remove line breaks spaces etc
	my $DNA_ic = anne_genomics::inverse_complement($DNA) ;
	read_condontable();
	my %results ;  # keys is the frame number 1 to 6
	for (my $i=0; $i<=2; $i++) {
		$results{$i+1} = translate_2_orfs(substr $DNA, $i) ;
		$results{$i+4} = translate_2_orfs(substr $DNA_ic, $i) ;
	}
	return %results ;
}

sub translate {
	my $DNA = shift ; 
	my $protein ;
	for ( my $i=0; $i<(length($DNA)-2); $i+=3) {
		my $codon=substr($DNA,$i,3);
		if (defined $codontable{$codon}) { $protein.=$codontable{$codon} }	else { $protein.='x' ;	}
	}
	return $protein ;
}


sub six_frames_orfs {
	# translate DNA to 6 frames
	my $DNA = uc(shift) ;
	$DNA =~ s/(\s|\n|\r|\x0d)//g; # remove line breaks spaces etc
	my $DNA_ic = anne_genomics::inverse_complement($DNA) ;
	read_condontable();
	my %results ;  # keys is the frame number 1 to 6
	for (my $i=0; $i<=2; $i++) {
		$results{$i+1} = translate_2_orfs(substr $DNA, $i) ;
		$results{$i+4} = translate_2_orfs(substr $DNA_ic, $i) ;
	}
	return %results ;
}

sub translate_2_orfs {
	my $DNA = shift ; 
	my $protein ;
	my $stopcodon = 1 ;
	for ( my $i=0; $i<(length($DNA)-2); $i+=3) {
		my $codon=substr($DNA,$i,3);
		if ($stopcodon and $codon =~ m/ATG|TTG|GTG/) { $stopcodon = 0 ; }
		if (defined $codontable{$codon} and !$stopcodon) { 
			$protein.=$codontable{$codon} 
		}	else { 
			$protein.='x' ;	
			$stopcodon = 1 ;
		}
	}
	return $protein ;
}



## the mandatory one (without it no package!!!)
1
