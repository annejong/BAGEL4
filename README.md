# BAGEL4
BAGEL4 v.1.2  June 2022

To download all dependencies use http://ngs.molgenrug.nl/bagel4/

---------- SOFTWARE dependencies ---------- 

perl5, <br>
hmmsearch3, <br>
blastall, <br>
glimmer3, <br>
librsvg2-bin,<br> 
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
																								<br>
---------- LIBRARIES ----------                                                                 <br>
																								<br>
Add the BAGEL4 library to PERL5LIB: e.g., export PERL5LIB=/data/pg-molgen/software/bagel4/lib   <br>
Or copy the files from the bagel4/lib to your perl5 lib folder ( e.g., /usr/lib/perl5)          <br>
																								<br>
---------- TOOLS ------------------                                                             <br>
																								<br>
Make symbolic link in MOODS folder in the BAGEL4 program dir                                    <br>
cd tools/MOODS-python-1.9.3/scripts                                                             <br>
ln -s ../MOODS/                                                                                 <br>
																								<br>
																								<br>
---------- UPDATING BAGEL ----------                                                            <br>
																								<br>
https://github.com/annejong/BAGEL4                                                              <br>
																								<br>
---------- EXECUTING BAGEL ----------                                                           <br>
																								<br>
To run BAGEL:                                                                                   <br>

bagel4_wrapper.pl is a wrapper around all files in a folder or subfolder and also will screens packed files (zip, rar) or folders recursively based on a regular expression
																																									   <br>
Without parameters bagel4_wrapper.pl the options will be listed:                                                                                                       <br>
			-s Sessiondir [default=./my_session]                                                                                                                       <br>
			-i queryfolder or zip_file containing DNA files in the FASTA format                                                                                        <br>
			-webserver [default=0; Meaning that the analysis will de done on the folder given after -i, else a query from the webserver is expected]                   <br>
			-glimmer make glimmer model [default=1; means that a global model will be made] Usefull for multiple entry FastA files with relative small DNA fragments   <br>
			-r regular expression for identifing files [default=.* meaning all files ]                                                                                 <br>
																																									   <br>
																																									   <br>


---------- EXAMPLES  ----------                                                                   <br>
																								  <br>
most basic option to execute BAGEL                                                                <br>
																								  <br>
e.g.  ./bagel4_wrapper.pl -i /var/genomes/Lactococcus_lactis_cremoris_MG1363_uid58837 -r \.fna\$  <br>
																								  <br>
e.g.  ./bagel4_wrapper.pl -i /var/dna/myfolder -r \.fna\$                                         <br>
																								  <br>
																								  <br>
================================                                                                  <br>
																								  <br>
Anne de Jong and Auke van Heel                                                                    <br>
university of Groningen                                                                           <br>
Nijenborgh 7                                                                                      <br>
9743 AG  Groningen                                                                                <br>
the Netherlands                                                                                   <br>
anne.de.jong@rug.nl                                                                               <br>
																								  <br>
				June 2022                                                                         <br>
																								  <br>
================================                                                                  <br>
