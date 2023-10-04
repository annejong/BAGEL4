package anne_misc ;

##############################################################

## 	Anne de Jong
##	University of Groningen
##	the Netherlands
##	anne.de.jong@rug.nl

##############################################################
##
##   
##  2013 oktober,  misc module 
##

use strict ;
use warnings ;

sub get_fonts {
	my %fonts ;
	my $folder = "/usr/molgentools/fonts";
	$fonts{FreeSerif} 			= "$folder/FreeSans.ttf";	
	$fonts{FreeSerifBold} 		= "$folder/FreeSerifBold.ttf";	
	$fonts{FreeSerifBoldItalic} = "$folder/FreeSerifBoldItalic.ttf";	
	$fonts{FreeSerifItalic} 	= "$folder/FreeSerifItalic.ttf";	
	$fonts{UbuntuRegular} 		= "$folder/Ubuntu-R.ttf";	
	$fonts{UbuntuRegularItalics}= "$folder/Ubuntu-RI.ttf";
	# below the default fonts
	$fonts{normal} 				= "$folder/Ubuntu-R.ttf";	
	$fonts{italics}				= "$folder/Ubuntu-RI.ttf";	
	$fonts{bold}				= "$folder/Ubuntu-B.ttf";	
	$fonts{bolditalics}			= "$folder/Ubuntu-BI.ttf";	
	return %fonts ;
}	

sub str_in_array {
	my ($str, @array) = @_ ;
	return grep { $_ eq $str } @array ;
}

sub repeat {
	my ($count, $str) = @_ ;
	print "===> $count, $str\n";
	my $result = '';
	for (my $i=1; $i<=$count; $i++) { $result .= $str ; }
	return $result ;
}

sub trim {
	 my $str = shift ;
	 $str =~ s/^\s+|\s+$//g ;
	 return $str;
}

sub de_space_array {
	# remove tailing spaces
	# replace internal space by _
	my @items = @_ ;
	for (my $i=0; $i< scalar(@items); $i++) {
		$items[$i] =~ s/^\s+|\s+$//g;
		$items[$i] =~ s/\ /\_/g ;
	}
	return @items ;
}


sub deg2rad {
	# convert degree values to radian
	my $deg = shift ;
	my $pi = 3.14159265358979 ;
	my $rad = ( $deg/360 ) * 2 * $pi ;
	return $rad ;
}

sub is_realnumber {
	my $number = shift ;
	my $result = 0 ;
	#$result = 1 if ($number =~ m/^-?\d+\.?\d*$/) ;# is a real number
	$result = 1 if ($number =~ m/^[-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+(?:\.\d*)?)?$/) ;# is a real number or scientific
	return $result ;
}

sub is_numeric {
	my $number = shift ;
	my $result = 0 ;
	$result = 1 if ($number =~ m/\d+/) ; 
	return $result ;
}

sub is_numeric_array {
	# check if all items are numeric
	my $result = 1 ;
	foreach my $number ( @_ ) {	if (!is_realnumber($number)) { $result = 0 ; last; } }
	return $result ;
}

sub log_array {
	my @result ;
	foreach my $number ( @_ ) {	push @result, log($number) ; }
	return @result ;

}
sub ceil {
	int($#_ + 0.99);
}

sub sum {
	my $sum =0 ;
	$sum += $_ for(@_);
	return $sum ;
}

sub median { $_[0]->[ @{$_[0]} / 2 ] }

sub median2
{
    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    if($len%2) #odd?
    {
        return $vals[int($len/2)];
    }
    else #even
    {
        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
}

sub average {
	my @numbers = @_ ;
	my $total = 0 ;
    foreach my $number (@numbers) { $total += $number; }
    my $average = $total / (scalar @numbers) ;
}	

sub log_N
{
  my $num = shift;
  my $base = shift;
  return log($num)/log($base);
}

sub roundup {
    my $n = shift;
    return(($n == int($n)) ? $n : int($n + 1))
}

sub count_str_in_str {
	# count the number of regex in a string
	my $regex = shift ;
	my $string = shift ;
	return () = $string =~ /$regex/g;
}
	


## the mandatory one (without it no package!!!)
1