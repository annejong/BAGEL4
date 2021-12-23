BAGEL4 v.1.0  July 2018

---------- SOFTWARE dependencies ---------- 

perl
hmmsearch3
blastall
glimmer3
librsvg2-bin
PfamScan     # this program is included in the BAGEL4.tar or download the latest version from  ftp://ftp.sanger.ac.uk/pub/databases/Pfam/Tools/


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

unpack BAGEL4.tar in e.g. /usr/bagel4 


---------- CONFIGURATION ---------- 

Open the bagel4.conf 
- change the PFAM database location if needed
- change the default number of cpu's (threads)
- for optimal performance of BAGEL keep the other setting untouched

LIBRARIES
- if BAGEL is installed in an other folder then /usr/bagel4 do the following:
	- open bagel4.pl to set proper path for the BAGEL4 program
	- edit line: use lib "/usr/bagel4/lib";
	- edit line: my %conf = read_conf('/usr/bagel4/bagel4.conf') ;


---------- TOOLS ------------------

Make symbolic link in MOODS folder
cd tools/MOODS-python-1.9.3/scripts
ln -s ../MOODS/ 


---------- UPDATING BAGEL ----------


---------- EXECUTING BAGEL ----------

There are two program to run BAGEL:

1. bagel4.pl is the main program and can handle single files

		bagel4.pl
			-s Sessiondir [default=current folder]
			-i input file should be a DNA file in the FASTA format
			-r NCBI reference genome
			-is_circular [default=1 for yes]
			-glimmer [default=0, 1 means that the training model is present in the sessiondir; normally made by bagel3_wrapper.pl]
			-graphics [default=1, 1 means that a graphic is drawn of the AOI]
			-o output file [default=results.txt]

2. bagel4_wrapper.pl is a wrapper around bagel4.pl and screens packed files (zip, rar) or folders recursively based on a regular expression


		bagel4_wrapper.pl
			-s Sessiondir [default=./my_session]
			-i queryfolder or zip_file containing DNA files in the FASTA format 
			-webserver [default=0; Meaning that the analysis will de done on the folder given after -i, else a query from the webserver is expected]
			-glimmer make glimmer model [default=1; means that a global model will be made] Usefull for multiple entry FastA files with relative small DNA fragments
			-r regular expression for identifing files [default=.* meaning all files ]


NOTE: executing bagel4.pl or bagel4_wrapper.pl without parameters will give you the list of options needed to run these programs



---------- EXAMPLES  ----------

most basic option to execute BAGEL

e.g.  ./bagel4.pl -i /var/genomes/Lactococcus_lactis_cremoris_MG1363_uid58837/NC_009004.fna 
e.g.  ./bagel4_wrapper.pl -i /var/genomes/Lactococcus_lactis_cremoris_MG1363_uid58837 -r \.fna\$

e.g.  ./bagel4_wrapper.pl -i /var/dna/myfolder -r \.fna\$ 


================================

Anne de Jong and Auke van Heel
university of Groningen
Nijenborgh 7
9743 AG  Groningen
the Netherlands
anne.de.jong@rug.nl

				July 2018

================================
