cd /data/bagel4/db_hmm/Bacteriocin

function my_hmmbuild {
	hmmbuild -n $1 $1.hmm $1.fasta
}		

my_hmmbuild bacteriocinII
my_hmmbuild CyaG
my_hmmbuild ggmotif
my_hmmbuild httB
my_hmmbuild labKC
my_hmmbuild lanl
my_hmmbuild leaderLanBC
my_hmmbuild leaderLanM
my_hmmbuild LE-DUF
my_hmmbuild LE-LAC481
my_hmmbuild LE-LanBC
my_hmmbuild LE-MER+2PEP
my_hmmbuild LinL
my_hmmbuild MA-2PEPA
my_hmmbuild MA-2PEPB
my_hmmbuild MA-DUF
my_hmmbuild MA-EPI
my_hmmbuild MA-LAC481
my_hmmbuild MA-NIS+EPI
my_hmmbuild MA-NIS
my_hmmbuild Sac_1
my_hmmbuild lasso


# concatenate all hmm models to one file <== depreciated

cat Bacteriocin/*.hmm > bagel4_bacteriocin_hmms
cat PF01721.hmm >>bagel4_bacteriocin_hmms
cat PF02052.hmm >>bagel4_bacteriocin_hmms
cat PF03047.hmm >>bagel4_bacteriocin_hmms
cat PF04369.hmm >>bagel4_bacteriocin_hmms
cat PF04454.hmm >>bagel4_bacteriocin_hmms
cat PF04604.hmm >>bagel4_bacteriocin_hmms
cat PF08109.hmm >>bagel4_bacteriocin_hmms
cat PF08129.hmm >>bagel4_bacteriocin_hmms
cat PF08130.hmm >>bagel4_bacteriocin_hmms
cat PF09221.hmm >>bagel4_bacteriocin_hmms
cat PF09683.hmm >>bagel4_bacteriocin_hmms
cat PF10439.hmm >>bagel4_bacteriocin_hmms
cat PF11758.hmm >>bagel4_bacteriocin_hmms
cat PF16934.hmm >>bagel4_bacteriocin_hmms


cd /data/bagel4/tables
# Add the motif(s) to bacteriocin_hmm.txt 

