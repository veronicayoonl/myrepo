# myrepo
making a change from Rstudio on my computer.
I am making this change from Github in the browser.
#Final project description CPSC565
 
#final_project_description.pl

#starting material -> 9hk_rawdata.fastq & 9rb_rawdata.fastq

#Objectives:
1) trim illumina universal adapters with Trimmomatic software
	command line: java -jar <path to trimmomatic> SE -threads 5 -phred33 -trimlog \n
	9hk_CGATGT.trimlog.txt /home/a-m/cpsc565_stud23/raw-data/9hk_CGATGT_L005_R1_001.fastq \n
	9hk_CGATGT.trimmed.fastq ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 \n
	LEADING:28 TRAILING:28 MINLEN:36
2) check how many lines and characters there are in each file (line_count.pl)
3) extract lines of nucleotide sequences with next if $line =~ /[@B+]/ \n
	(get rid of Illumina Phred scores) fastq_to_fasta.pl
4) discard 6 first nucleotides of reads (because of high kmer content based upon FASTQC \n
	report) minus_6.pl
5.1) re-check number of lines and characters (for comparison with raw-data)
5.2) Measure GC content in both files (gc_test.pl): this is done just to see if GC \n
	composition matches that of the FASTQC report.
6) make line lengths divisible by 3 (go through every line, if any line length is not \n
	divisible by three concatenate with next line until it is divisible by three) \n
	chomp_seq_alt.pl
7) translate three-based nucleotides to amino-acids (stop codons are represented by "*" \n
	symbol) with module Gencode.pm
7.1) Execute package Gencode.pm with gencode.pl (warning: this takes more than an hour)
8) split protein sequence by * and then incorporate new line character (\n) prot_ext.pl
9) discard protein sequences that are composed of less than 30 amino acids (threshold of \n	
	30 was arbitrarily defined)
10)incorporate header to each file (with > Protein_counter) --> this is to convert text \n	
	to fasta file (9 and 10 are executed by prot_align.pl)
11) because output file name lengths progressively increased, in the end these were \n
	name-modified (For example: mv 9hk*.prot_align 9hk.final.fasta)

 check against pfam database for hmm hits:
 1) module load hmmer
 2) database is in path /home/classroom/cpsc565/Pfam-A.hmm
 3) command line: hmmsearch /home/classroom/cpsc565/Pfam-A.hmm 9hk.final.fasta > 9hk.out \n
 	(warning: this takes more than an hour!)
 4) command line: hmmsearch /home/classroom/cpsc565/Pfam-A.hmm 9rb.final.fasta > 9rb.out\n	
	(warning: this takes more than an hour!)


CODES:

TO COUNT LINES AND CHARACTERS

#!/usr/bin/perl
#line_count.pl
use strict; use warnings;

die "Usage: line_count.pl <file>\n" if (@ARGV != 1);
print "Processing file $ARGV[0]\n";

my $lines = 0;  #initialize count to zero
my $characters = 0;

while (my $line = <>) {
        $lines++;
        $characters += length($line);  #this is equal to $character = \n
        							   #$character + length($line), this will make \n
        							   #$characters to count all characters of the input \n
        							   #file 
}
print "$lines\t$characters\n"; #prints number of lines tab number of characters


FASTQ TO "FASTA" 

#!/usr/bin/perl
# fastq_to_fasta.pl
use strict; use warnings;

my ($infile) = @ARGV;

my $outfile = "$infile.fasta1";

open(my $in, "<$infile") or die "error reading $infile. $!";
open(my $out, ">$outfile") or die "error creating $outfile.$!";

while (my $line = <$in>) {
        next if ($line =~ m/[B@+]/gi);# proceed to next line if you find B, @, or + \n
        							  #symbols at each line
        print $out "$line";    #print $line in output file
        print length($line);
        print "\n";
        my $line++;   #line is incremented by one during each loop iteration
}

close $in;
close $out;



REMOVE FIRST 6 NUCLEOTIDES

#!/usr/bin/perl
# fastq_to_fasta.pl
use strict; use warnings;

my ($infile) = @ARGV;

my $outfile = "$infile.m.3.fasta";


open(my $in, "<$infile") or die "error reading $infile. $!";
open(my $out, ">$outfile") or die "error creating $outfile.$!";

while (my $line = <$in>) {
        substr($line, 0, 6)="";  #take six first nucleotides/characters and convert to "".
        print $out "$line";
}

close $in;
close $out;


CALCULATE GC COMPOSITION 

#!/usr/bin/perl
#gc_test.pl
use strict; use warnings;


my ($infile) = @ARGV;
open(my $in, "<$infile") or die "error reading $infile. $!";

my $G=0;
my $H=0;
my $C=0;
my $Ratio=0;
while (my $line = <$in>) {
        $G=$G+$line =~ tr/G/G/;   #also could have used $G += $line
        $C=$C+$line =~ tr/C/C/;
        $H=$H+length($line);
        $Ratio =($G+$C)/$H;

}

print "GC composition is $Ratio \n";

close $in;


MAKE LINE LENGTH DIVISIBLE BY 3

#!/usr/bin/perl
# chomp_seq.pl
use strict; use warnings;

my ($infile) = @ARGV;

my $outfile = "$infile.csh";


open(my $in, "<$infile") or die "error reading $infile. $!";
open(my $out, ">$outfile") or die "error creating $outfile.$!";

my $asdf = "";
while (<$in>) {
        chomp;
		$asdf = ($_).$asdf; #concatenate $asdf with new line 
        if (length($asdf)%3 == 0) {    #if concatenation of lines is divisible by three 
        print $out "$asdf\n";   #print line 
        $asdf = ""; #reset concatenation of lines to ""
        }
}
print $out "$asdf\n"; #print last concatenation

close $in;
close $out;


TRANSLATE CODONS TO AMINO ACIDS

package Gencode;
use strict; use warnings;
my %GENCODE = (
        'AAA' => 'K', 'AAC' => 'N', 'AAG' => 'K', 'AAT' => 'N',
        'ACA' => 'T', 'ACC' => 'T', 'ACG' => 'T', 'ACT' => 'T',
        'AGA' => 'R', 'AGC' => 'S', 'AGG' => 'R', 'AGT' => 'S',
        'ATA' => 'I', 'ATC' => 'I', 'ATG' => 'M', 'ATT' => 'I',
        'CAA' => 'Q', 'CAC' => 'H', 'CAG' => 'Q', 'CAT' => 'H',
        'CCA' => 'P', 'CCC' => 'P', 'CCG' => 'P', 'CCT' => 'P',
        'CGA' => 'R', 'CGC' => 'R', 'CGG' => 'R', 'CGT' => 'R',
        'CTA' => 'L', 'CTC' => 'L', 'CTG' => 'L', 'CTT' => 'L',
        'GAA' => 'E', 'GAC' => 'D', 'GAG' => 'E', 'GAT' => 'D',
        'GCA' => 'A', 'GCC' => 'A', 'GCG' => 'A', 'GCT' => 'A',
        'GGA' => 'G', 'GGC' => 'G', 'GGG' => 'G', 'GGT' => 'G',
        'GTA' => 'V', 'GTC' => 'V', 'GTG' => 'V', 'GTT' => 'V',
        'TAA' => '*', 'TAC' => 'Y', 'TAG' => '*', 'TAT' => 'Y',
        'TCA' => 'S', 'TCC' => 'S', 'TCG' => 'S', 'TCT' => 'S',
        'TGA' => '*', 'TGC' => 'C', 'TGG' => 'W', 'TGT' => 'C',
        'TTA' => 'L', 'TTC' => 'F', 'TTG' => 'L', 'TTT' => 'F'
);
sub translate_codon {     #this function returns a codon that is passed to the next \n
						  #function to be used as a look-up key in the Gencode hash
        my ($codon) = @_;
        if (exists $GENCODE{$codon}) {return $GENCODE{$codon}}
        else                         {return 'X'}
}
sub translate_sequence {
        my ($seq) = @_;
        my $pep;
        #my $codon;
        my $k = length($seq);
        for (my $i = 0; $i < length($seq); $i += 3) {#for loops have three parameters;\n
        										  #initialization, validation, and update
                my $codon = substr($seq, $i, 3);
                print "$codon\n"; 
                print "$i/$k\n";  
        if (exists $GENCODE{$codon}) {$pep .= $GENCODE{$codon}}
        else                         {$pep .= 'X'}
        }
return $pep;
}
1;


to execute Gencode.pm module use:

#!/usr/bin/perl
#gencode.pl
use strict; use warnings;
use Gencode;

my ($infile) = @ARGV;

my $outfile = "$infile.prot";

open(my $in, "<$infile") or die "error reading $infile. $!";
open(my $out, ">$outfile") or die "error creating $outfile.$!";
my $k;

while (my $line = <$in>) {
        print length($line);
        print "\n";
		chomp;
		my $pep = Gencode::translate_sequence($line);
        print $out $pep;
}

close $in;
close $out;




EXTRACT PROTEIN SEQUENCES THAT ARE BETWEEN ASTERISKS (STOP CODONS)

#!/usr/bin/perl
# chomp_seq.pl
use strict; use warnings;

my ($infile) = @ARGV;

my $outfile = "$infile.protext";


open(my $in, "<$infile") or die "error reading $infile. $!";
open(my $out, ">$outfile") or die "error creating $outfile.$!";

while (<$in>) {
        my @prots = split /\Q*\E/, $_;#\Q and \E are used to ignore * as special character
        for (@prots) {
                print $out "$_\n";
        }
}

close $in;
close $out;

REMOVE PROTEINS THAT ARE COMPOSED OF LESS THAN 30 AMINOACIDS AND THEN INCORPORATE HEADER\n
(>PROTEIN_COUNTER)

#!/usr/bin/perl
#prot_align.pl
use strict; use warnings;

my ($infile) = @ARGV;

my $outfile = "$infile.ready_align";

open(my $in, "<$infile") or die "error reading $infile. $!";
open(my $out, ">$outfile") or die "error creating $outfile.$!";

my $counter = 1; #counter begins at 1
while (my $line = <$in>) {
        if (length($line) >= 30) {  #only print those lines that are above or equal to 30
        print $out ">protein_$counter\n";  #print header ">protein_counter" new line \n
        								   #character 
        print $out "$line";
        $counter = $counter + 1; #increment counter by 1
        }
}

close $in;
close $out;


