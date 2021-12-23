
function update_pfam {
	wget http://pfam.xfam.org/family/$1/hmm -O $1.hmm
}

# added for BAGEL4
update_pfam PF16934
update_pfam PF09221
update_pfam PF10439
update_pfam PF11758
update_pfam PF03047

# inherited from bagel3
update_pfam PF00005
update_pfam PF00069
update_pfam PF00072
update_pfam PF00082
update_pfam PF00296
update_pfam PF00486
update_pfam PF00512
update_pfam PF00535
update_pfam PF00733
update_pfam PF00881
update_pfam PF01320
update_pfam PF01721
update_pfam PF01944
update_pfam PF02052
update_pfam PF02441
update_pfam PF02518
update_pfam PF02624
update_pfam PF02794
update_pfam PF03412
update_pfam PF03526
update_pfam PF03857
update_pfam PF04055
update_pfam PF04369
update_pfam PF04454
update_pfam PF04604
update_pfam PF04738
update_pfam PF05147
update_pfam PF06580
update_pfam PF08109
update_pfam PF08129
update_pfam PF08130
update_pfam PF08951
update_pfam PF09221
update_pfam PF09683
update_pfam PF10439
update_pfam PF12730
update_pfam PF13437
update_pfam PF13471
update_pfam PF13575
update_pfam PF14028
update_pfam PF04738


