package bagel4_graphics ;

#     drawing module for BAGEL  
# 	
#	University of Groningen
#	the Netherlands
#	anne.de.jong@rug.nl
#
#   
#  2012 December


use strict;
use warnings;
use Storable ;
use GD;
use GD::Polyline;

BEGIN {
    use Exporter ();
}

my $width = 11000;					
my $height = 3000;	
my $canvasmargin = 10;
my $topmargin = 1500 ;
my $linespacing = 50 ; 
my $genesize = 400 ;
my $bpline = 20000;					# Number of bp per line
my $bpixel = $bpline / $width ;		# Number of bp per pixel
my $font = "/usr/share/fonts/truetype/freefont/FreeSans.ttf";
my $fontsize = 90 ;
my $image ;
my %kleur ;
my $prev_genename_xpos ;


sub draw_AOI {
	my ($filename, $infont, $queryLen, %table) = @_ ;
	$font = $infont ;
	$image = new GD::Image ( $width+$canvasmargin,$height+$canvasmargin );
	# http://www.rapidtables.com/web/color/RGB_Color.htm
	$kleur{white} 	= $image->colorAllocate( 255, 255, 255 );
	$kleur{red} 	= $image->colorAllocate( 255, 0, 0 );
	$kleur{greenI} 	= $image->colorAllocate( 0, 204, 0 );
	$kleur{greenII} = $image->colorAllocate( 0, 153, 0 );
	$kleur{greenIII}= $image->colorAllocate( 0, 102, 0 );
	$kleur{purple} 	= $image->colorAllocate( 76, 0, 153 );
	$kleur{pink} 	= $image->colorAllocate( 255, 153, 255 );
	$kleur{greenblue}= $image->colorAllocate( 0, 153, 153 );
	$kleur{lime} 	= $image->colorAllocate( 0, 255, 0 );
	$kleur{blue} 	= $image->colorAllocate( 0, 0, 255 );
	$kleur{yellow} 	= $image->colorAllocate( 204, 204, 0 );
	$kleur{black} 	= $image->colorAllocate( 0, 0, 0 );
	$kleur{gray} 	= $image->colorAllocate( 192, 192, 192 );
	ruler($queryLen);
	$prev_genename_xpos = 0 ;
	foreach my $key (sort {$table{$a}{Gene_start} <=> $table{$b}{Gene_start}} keys %table) {
		my $color = $table{$key}{color} ;
		my $genename = $key.' '.$table{$key}{pfamname} ;
		$genename = $key.'_'.$table{$key}{contextname} if ($table{$key}{contextname} ne '' ) ;
		draw_one_gene($genename,$table{$key}{Gene_start},$table{$key}{Gene_end},$table{$key}{strand}, $kleur{$color}, $table{$key}{length}) ;
	}
	save_image($filename) ;
}



sub ruler {
	my $queryLen = shift ;
	$image->setThickness(3) ;
	my $x2 = $canvasmargin + ($queryLen/$bpline) * ($width- (2*$canvasmargin) ) ;
	my $y1 = $topmargin + $genesize/2 ;
	# $x1,$y1,$x2,$y2,$color
	$image->line($canvasmargin, $y1 , $x2, $y1 ,$kleur{black});
	$image->setThickness(4) ;
	# ruler
	$y1 = 0.8 * $height ; 
	my $y2 = $y1 - 50 ;
	my $angle = 0 ;
	$image->line($canvasmargin,$y1, $x2,$y1 ,$kleur{black});
	for ( my $count = 0; $count <= 18000; $count = $count+2000) {
		#print "Count=$count\n";
		my $x1 = bp2_xpos($count) + $canvasmargin ;
		$image->line($x1, $y1, $x1,$y2 ,$kleur{black});
		$image->stringFT($kleur{black},$font,$fontsize,$angle, $x1 , $y2, $count) ;
	}
	$image->setThickness(1) ;
}

sub draw_one_gene {
# calculation the gene positions on the canvas
    my ($genename, $start, $end, $strand, $color, $genelen) = @_ ;
	my $x 	= $start / $bpline ;  # calculate to which line the gene should be drawn  
    $x =~ s/.*\.//;	# get the decimals
    $x = "0.$x";
	$x = $x * $width  ;
    my $y 	= $topmargin + $linespacing * (int ($start / $bpline )); 
	my $l	= ($end-$start) / $bpixel ;

	# print the gene name
	my $angle = 0.62 ;
	my $x_shift = 5 ;
	$x_shift = 100 if ($genelen > 100) ;
	$x_shift = 200 if ($genelen > 200) ;
	$x_shift = 400 if ($genelen > 800) ;
	if (!($genelen < 75 and $color eq '12')) { 				# do not write gene names for small gray genes to prevent overlapping text $genename
		my $genename_xpos = $x + $x_shift ;
		$prev_genename_xpos = $prev_genename_xpos + 200 ; 	# this is added to prevent overlapping gene names
		$genename_xpos = $prev_genename_xpos if ($genename_xpos < $prev_genename_xpos ) ;
		$prev_genename_xpos = $genename_xpos ;
		if ($genename =~ m/.*\;(.*)/g) {
			$image->stringFT($kleur{black},$font,$fontsize,$angle,$genename_xpos, $y-16, $1) 
		}
	}	
	# draw the gene and check for end of the canvas
	if ($x+$l > $width) {
		if ($strand eq '+') { # this is the upper strand 1
			$image->filledRectangle ($x,$y,$width,$y+10, $color) ;
			draw_gene_polygon (1, $y+$linespacing, $x+$l-$width, $genelen, $strand, $color) ;
		} else { # this is the lower strand -1
			$image->filledRectangle (1,$y+$linespacing,($x+$l)-$width,$y+$linespacing+10, $color) ;
			draw_gene_polygon ($x, $y, $width-$x, $genelen, $strand, $color) ;
		}
	} else {  # just draw the complete gene
		draw_gene_polygon ($x, $y, $l , $genelen, $strand, $color) ;
    }
	#print "Coords in pixels:\t".$x."\t".$y."\t".$l."\n";
}

sub draw_gene_polygon {
# drawing the gene polygon
	my ($x, $y, $l, $genelen, $strand, $color) = @_ ;
    my $pointsize = $genesize / 2 ;
	$pointsize = $genesize / 4 if ($l < 200) ;
	$pointsize = $genesize / 7 if ($l < 100) ;
	$pointsize = $genesize / 10 if ($l < 70) ;
	my $polygon = new GD::Polyline;
	if ($strand eq '+' ) {
		$polygon->addPt( $x,$y);
		$polygon->addPt( $x+$l-$pointsize,$y);
		$polygon->addPt( $x+$l, $y+$genesize/2);
		$polygon->addPt( $x+$l-$pointsize,$y+$genesize);
		$polygon->addPt( $x,$y+$genesize );
    } elsif ($strand eq '-' ) {
		$polygon->addPt( $x+$l,$y);
		$polygon->addPt( $x+$pointsize,$y);
		$polygon->addPt( $x, $y+$genesize/2);
		$polygon->addPt( $x+$pointsize,$y+$genesize);
		$polygon->addPt( $x+$l,$y+$genesize );
    }
	$image->filledPolygon( $polygon, $color );
}	

sub bp2_xpos {
# calculated the x pixel postion on basis of the bp position
    my $start = shift ;
	my $x 	= $start / $bpline ;  # calculate the which line the gene should be drawn  
    $x =~ s/.*\.//;	# get the decimals
    $x = "0.$x";
	$x = $x * $width  ;
	return $x ;
}
 
sub bp2_ypos { 
# return the y position on basis of bp position
	my $start = shift ;
	my $y 	= $topmargin + $linespacing * (int ($start / $bpline )); 
	return $y ;
}

    
sub save_image {
	my $file = shift ;
	open( PICT, ">$file" ) or die( "Cannot save $file" );
	binmode( PICT );
	print( PICT $image->png() );
	close PICT ;
}

1;