# Build the hmm models for the promoter prediction

hmmbuild  --dna -n Sigma_N16 SigmaA_N16.hmm SigmaA_N16.aln
hmmbuild  --dna -n Sigma_N17 SigmaA_N17.hmm SigmaA_N17.aln
hmmbuild  --dna -n Sigma_N18 SigmaA_N18.hmm SigmaA_N18.aln
hmmbuild  --dna -n min35 min35.hmm min35.aln

cat SigmaA_N16.hmm SigmaA_N17.hmm SigmaA_N18.hmm min35.hmm > SigmaA.hmm
