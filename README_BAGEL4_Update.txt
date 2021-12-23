Update of Bagel4

1) Goto folder "db_proteins"
	update content of the files bagel4_class1bacteriocin_db.txt, bagel4_class2bacteriocin_db.txt and bagel4_class3bacteriocin_db.txt
	run: ./UPDATE_bagel4_bacteriocin_db.pl
	
2) After updating the db_proteins, goto folder "db_uniref90_xml"
	run: ./00.update_uniref_xml.pl
	
3) Goto folder "swissprot"
	run: 00.README.update_swissprot.sh

