
**mRNA_Markup**


1. Description
2. Installation requirements
3. Formatting a local nucleotide and peptide database using BLAST+
4. Integration with a local version of GALAXY

<br>
<br>
<ol>
<li>Description</li>
<br>
This Galaxy tool is based on BioExtract mRNA Markup workflow.
This workflow represents comprehensive annotation 
and primary analysis of a set of transcripts.

It consists of several "steps" each using several analytic tools:<br>
<br>
0. Submit input: initial mRNA (query file) <br>
		 bacterial contamination database in FASTA format (a file representing typical bacterial hosts) <br>
		 reference protein database in FASTA format (proteins most likely to have homologs in the mRNA translations of the input)<br>
		 comprehensive protein database in FASTA format (to be searched when the reference protein set did not give hits)<br>
1. Eliminate Vector Contamination<br>
2. Eliminate bacterial contamination<br>
3. Find matches in a reference protein database <br>
3.1. Identify potential full-length coding sequences <br>
3.2. Identify potential chimeric sequences<br>
4. Find matches in a Comprehensive Protein database<br>
5. Find matches in Protein Domain Database<br>
6. Produce summary report<br>
<br>
<br>
For more information see:<br>
(1) http://www.bioextract.org/scripts/index.jsp<br>
(2) http://bioservices.usd.edu/mrnamarkup.html<br>

<br>
<br>

<li>Installation requirements</li>
<br>
Before you start make sure that <br>
(1) ncbi-blast+<br>
(2) biopython<br>
(3) MuSeqBox<br>
are already installed on your computer.<br>
<br>
Note: Version 4.1 of MuSeqBox (available from 
http://www.plantgdb.org/MuSeqBox/MuSeqBox.php)
is based on BLAST+ applications and was tested with ncbi-blast-2.2.24, 
available from http://blast.ncbi.nlm.nih.gov/.  
As there are numerous changes in the output
format between BLAST+ and the legacy BLAST applications, MuSeqBox will fail
to parse output from older BLAST versions.
<br>


<br>
<li>Formatting a local nucleotide and peptide database using BLAST+</li>
<br>
The workflow uses a sample input file consisting of Arabidopsis mRNA and searches the following BLAST databases by default:<br><br>

BacteriaDB and RefProtDB	 	http://www.bioextract.org/Download?action=zip&file=AllOutput.zip&src=/usr/local/BioStreamServer/tmpFiles/mRNA_Markup_Start/1291915878836<br>
AllProtDB 				http://www.uniprot.org/uniref/?query=identity:0.9+taxonomy:33090&format=*&compress=yes<br>
Smart_LE/Smart				ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/big_endian/<br>
UniVec					ftp://ftp.ncbi.nih.gov/pub/UniVec/<br>

<br>
Command for formatting nucleotide databases:<br>
makeblastdb -dbtype nucl -in DATABASE -parse_seqids<br>
<br>

Command for formatting for protein databases:<br>
makeblastdb -dbtype prot -in DATABASE -parse_seqids<br>

<br><br>
Note: Smart databases do not need formatting. This databaes is already formatted and packed and available from the website mentioned above.

<br><br>
<li>Integration with a local version of GALAXY</li>
<br>
Simply copy the galaxy_dist directory from the repository onto the GALAXY directory tree.<br>
The mRNA_Markup was tested only with a local GALAXY version.<br>
</ol>

<br>
********************
Katarzyna Wręczycka<br>
email: kw292555@students.mimuw.edu.pl<br>
18.01.2013

