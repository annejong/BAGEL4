# BAGEL4
BAGEL4 v.1.2  June 2022

---------- SOFTWARE dependencies ---------- 

perl5
hmmsearch3
blastall
glimmer3
librsvg2-bin
PfamScan     # Download the latest version from  ftp://ftp.sanger.ac.uk/pub/databases/Pfam/Tools/


---------- DATABASE dependencies ---------- 

PFAM-A database:
	- HOW TO INSTALL PFAM database
		- get the latest release from: ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/   ....../Pfam-A.hmm.gz
			wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/Pfam-A.hmm.gz
			wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/Pfam-A.fasta.gz
			wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/Pfam-A.seed

		- Unzip the files
			gunzip Pfam*.gz

		- create index binaries
			hmmpress Pfam-A.hmm
			
		- Index Pfam-A.fasta for BLAST searches
			formatdb -i Pfam-A.fasta -p T


optional download: uniref_db = uniref50.fasta


---------- INSTALLING BAGEL ---------- 

unpack BAGEL4.tar 


---------- CONFIGURATION ---------- 

Open the bagel4.conf 
- Adapt folder locations to your system
- Optionally change the default number of cpu's (threads)
- For optimal performance of BAGEL keep the other setting untouched

---------- LIBRARIES ----------

Add the BAGEL4 library to PERL5LIB: e.g., export PERL5LIB=/data/pg-molgen/software/bagel4/lib
Or copy the files from the bagel4/lib to your perl5 lib folder ( e.g., /usr/lib/perl5)

---------- TOOLS ------------------

Make symbolic link in MOODS folder in the BAGEL4 program dir
cd tools/MOODS-python-1.9.3/scripts
ln -s ../MOODS/ 


---------- UPDATING BAGEL ----------

https://github.com/annejong/BAGEL4

---------- EXECUTING BAGEL ----------

To run BAGEL:

bagel4_wrapper.pl is a wrapper around all files in a folder or subfolder and also will screens packed files (zip, rar) or folders recursively based on a regular expression

Without parameters bagel4_wrapper.pl the options will be listed:
			-s Sessiondir [default=./my_session]
			-i queryfolder or zip_file containing DNA files in the FASTA format 
			-webserver [default=0; Meaning that the analysis will de done on the folder given after -i, else a query from the webserver is expected]
			-glimmer make glimmer model [default=1; means that a global model will be made] Usefull for multiple entry FastA files with relative small DNA fragments
			-r regular expression for identifing files [default=.* meaning all files ]




---------- EXAMPLES  ----------

most basic option to execute BAGEL

e.g.  ./bagel4_wrapper.pl -i /var/genomes/Lactococcus_lactis_cremoris_MG1363_uid58837 -r \.fna\$

e.g.  ./bagel4_wrapper.pl -i /var/dna/myfolder -r \.fna\$ 


================================

Anne de Jong and Auke van Heel
university of Groningen
Nijenborgh 7
9743 AG  Groningen
the Netherlands
anne.de.jong@rug.nl

				June 2022

================================
