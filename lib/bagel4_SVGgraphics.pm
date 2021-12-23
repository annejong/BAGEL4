package bagel3_SVGgraphics ;

#   SVG drawing module for BAGEL  
# 	
#	University of Groningen
#	the Netherlands
#	anne.de.jong@rug.nl
#
#   
#  2013 Februari


use strict ;
use warnings ;


my $href_pfam = 'http://pfam.sanger.ac.uk/family' ;
my $scale 	= 1 ;
my $canvas_width = '1200px' ; # the image resolution
my $canvas_height = '500px' ;
my $viewboxX = 24000 ;	# the viewport in the html page
my $viewboxY = 6000 ;
my $marginL	= 200 ; 
my $marginT = 200  ;
my $genewidth = 700 ;
my $fontsize = 200 ;
my $font = 'verdana';

my %table ;
my %colors ;
my @output ;


sub draw_SVG {
	%table = @_ ;
	%colors = allocate_rgb_colors() ;
	push @output, "<svg width=\"$canvas_width\" height=\"$canvas_height\" viewBox=\"0 0 $viewboxX $viewboxY\" xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\">" ;
	my $DNApostion = $marginT + $scale * 300 ;
	foreach my $key (sort keys %table) {
		my $genename = $key.' '.$table{$key}{pfamname} ;
		$genename = $key.'|'.$table{$key}{contextname} if ($table{$key}{contextname} ne '' ) ;
		push @output,  draw_gene($key, scaleY($DNApostion), $table{$key}{color}) ; 
		push @output,  draw_genename($key, scaleY($DNApostion),$genename) ; 
	}
	push @output,  draw_ruler();
	push @output,  "</svg>" ;
	return @output ;
}

sub draw_ruler {
	my $rulerypos = 2000 ;
	my $ticksize = 130 ;
	# the line
	my $x2 = $marginL + 20000 ;
	my @result = "\n\t<line x1=\"$marginL\" y1=\"$rulerypos\" x2=\"$x2\" y2=\"$rulerypos\" stroke=\"$colors{black}\" stroke-width=\"4\"/>" ;
	# the ticks
	for (my $i=0; $i <= 20000; $i=$i+2000 ) {
		my $x = $marginL + $i ;
		my $y = $rulerypos - $ticksize - 20 ;
		my $y2 = $rulerypos - $ticksize ; 
		push @result, "\t<line x1=\"$x\" y1=\"$rulerypos\" x2=\"$x\" y2=\"$y2\" stroke=\"$colors{black}\" stroke-width=\"3\"/>" ;
		$x = $x + 50 ;
		my $angle = "rotate(-90 $x, $y)" ;  
		push @result, "\t<text x=\"$x\" y=\"$y\" fill=\"$colors{black}\" font-size=\"200\" transform=\"$angle\">$i</text>" ;

	}	
	return @result ;
}



sub draw_genename {
	my ($key, $y, $genename) = @_ ;
	my $x_shift = 5 ;
	$x_shift = 100 if ($table{$key}{length} > 100) ;
	$x_shift = 200 if ($table{$key}{length} > 200) ;
	$x_shift = 400 if ($table{$key}{length} > 800) ;
	$x_shift += $table{$key}{Gene_start} ;
	my $x = scaleX($x_shift) ;
	$y = $y - $genewidth - 100;
	my $angle = "rotate(-30 $x,$y)" ;
	my $text = "\t<text x=\"$x\" y=\"$y\" fill=\"$colors{black}\" font-size=\"$fontsize\" font-family=\"$font\" transform=\"$angle\">$genename</text>" ;
	return $text ;
}

sub draw_gene {
	my ($key, $y, $color) = @_ ;
	my $x = scaleX($table{$key}{Gene_start}) ;
	my $l = scaleX($table{$key}{length}) ;
	my $p = 0.7 * $genewidth ;  # size of the point
	$p =$l if ($l < $p ) ; 
	# define the coords of the gene in + and - strand
	my @points ;
	if ($table{$key}{strand} =~ m/.*\-,*/g) {  # - or -1 of strand-
		push @points, $x+$l.','.$y;
		push @points, $x+$l.','.($y-$genewidth);
		push @points, $x+$p.','.($y-$genewidth);
		push @points, $x   .','.($y-$genewidth/2);
		push @points, $x+$p.','.$y;
	} else {
		push @points, $x.','.$y;
		push @points, $x+$l-$p.','.$y;
		push @points, $x+$l.','.($y-$genewidth/2);
		push @points, $x+$l-$p.','.($y-$genewidth);
		push @points, $x.','.($y-$genewidth);
	}
	# draw the gene polygon
	my @polygon = "\t<a xlink:href=\"href_pfam/$table{$key}{pfamname} target=_blank\">" ;
	push @polygon, "\t\t<polygon points=\"".join (" ",@points)."\" style=\"fill:$colors{$color};stroke:black;stroke-width:1\"/>" ;
	push @polygon, "\t</a>" ;
	return @polygon ;
}

sub scaleX {
	my $x = shift ;
	my $result = $marginL + $scale * $x ;
	return $result ;
}

sub scaleY {
	my $y = shift ;
	my $result = $marginT + $scale * $y ;
	return $result ;
}

sub allocate_rgb_colors {
	# http://www.rapidtables.com/web/color/RGB_Color.htm
	my %kleur ;
	$kleur{white}		='rgb(255,255,255)';
	$kleur{red}			='rgb(255,0,0)';
	$kleur{greenI}		='rgb(0,204,0)';
	$kleur{greenII}		='rgb(0,153,0)';
	$kleur{greenIII}	='rgb(0,102,0)';
	$kleur{purple}		='rgb(76,0,153)';
	$kleur{pink}		='rgb(255,153,255)';
	$kleur{greenblue}	='rgb(0,153,153)';
	$kleur{lime}		='rgb(0,255,0)';
	$kleur{blue}		='rgb(0,0,255)';
	$kleur{yellow}		='rgb(204,204,0)';
	$kleur{black}		='rgb(0,0,0)';
	$kleur{gray}		='rgb(192,192,192)';
	return %kleur ;
}

1;