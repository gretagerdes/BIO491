# Genome-Wide Association Study of SARS-CoV-2 Variants #
Genome-wide association studies (GWAS) are used to associate genetic variations with particular phenotypic qualities. In this repository, we will provide information about carrying out a GWAS for transmissibility in 4 SARS-CoV-2 variants. 
## Background Information ##
In order to gain a better understanding of the alignment algorithms used in this GWAS, we will first discuss the Needleman-Wunsch method, a small-global alignment technique for biological sequences. The algorithm takes an input of two strings (DNA sequences) and constructs a matrix with one string horizontally across the top row and the other vertically down the first column. For this example, we will use the strings: GAATTC and GATTA. 

<i></i>  | G | A | A | T | T | C 
--- | --- | --- | --- |--- |--- |--- 
<i></i> | <i></i> | <i></i> | <i></i> | <i></i> | <i></i> | <i></i>
**G**  | <i></i> | <i></i> | <i></i> | <i></i> | <i></i> | <i></i>
**A**  | <i></i> | <i></i> | <i></i> | <i></i> | <i></i> | <i></i>
**T**  | <i></i> | <i></i> | <i></i> | <i></i> | <i></i> | <i></i>
**T**  | <i></i> | <i></i> | <i></i> | <i></i> | <i></i> | <i></i>
**A**  | <i></i> | <i></i> | <i></i> | <i></i> | <i></i> | <i></i>

The values of the matrix are filled in based on the following scoring system: 
* Match: The two letters at the current index are the same (score = + 2) 
* Mismatch: The two letters at the current index are different (score = -1) 
* Gap: One letter is aligned to a gap in the other string (score = - 2) 
The top-left cell is filled with a 0 and the scoring system is used to fill in the remaining cells. Once the matrix is filled, we get the following result: 

<i></i>  | G | A | A | T | T | C 
--- | --- | --- | --- |--- |--- |--- 
0 | -2 | -4 | -6 | -8  | -10 | -12
**G**  | -2 | 2 | 0 |  -2 | -4 | -6 | -8
**A**  | -4 | 0 | 4 | 2 | 0 | -2 | 4
**T**  | -6 | -2 | 2 | 3 | 4 | 2 | 0
**T**  | -8  | -4| 0 | 1 | 5 | 6 | 4
**A**  | -10 | -6| -2 | 2 | 3 | 4 | 5

We created a MATLAB function that takes in an array of the two sequence strings, the completed scoring matrix, and the scoring values, and outputs the correct alignment of the sequences. The traceback algorithm works backwards through the values of the matrix to determine the alignment. We also included the MATLAB workbook that contains the array of the two sequences, completed scoring matrix, and scoring values for this example. Those files can be found in the MATLAB folder of this repository. 
For this examples, the correct alignment of the sequences is: 
G A A T T C and G – A A T A. 

## GWAS Pipeline ##
For the next section of this repository, we will layout the GWAS pipeline for the SARS-CoV-2 data. The genetic data from the SARS-CoV-2 wild type and four variants were collected from the National Center for Biotechnology Information virus SARS-CoV-2 Data Hub. The FASTA files used in this study are located in the “DNA files” folder in this repository. The following code can be executed in Ubuntu. 
The reference genome used in this study is GCF_009858895.2_ASM985889v3_genomic.fna. 
The sequence reads from the variant genomes are created using ART, a series of simulation tools that generate synthetic, illumina sequencing reads: 
```bash
sudo apt-get install art_illumina
```
Documentation for [art illumina](https://manpages.debian.org/stretch/art-nextgen-simulation-tools/art_illumina.1). 
For this GWAS, we will be creating paired end reads ( -p) with read lengths of 150 base pairs, 30 fold coverage, a mean size of DNA fragments of 200 base pairs ( -m 200) , and a standard deviation of 10 base pairs ( -s 10). We will also be using HS25 Illumina sequencing system ( -ss HS25). 
```bash
art_illumina -i (variant FASTA) -p -l 150 -ss HS25 -f 30 -m 200 -s 10 -o paired_data_(variant name) 
```
Next, we will use the [Burrows-Wheeler Aligner](http://manpages.ubuntu.com/manpages/bionic/man1/bwa.1.html) (BWA) software package was used to build an index for the reference genome. 

```bash
sudo apt-get install bwa
```
Index the reference with headers called “Ref”: 
```bash
bwa index -p Ref
```
Use BWA to align the variant illumina reads against the indexed reference genome. Complete the command twice for the paired-end reads: 
```bash
bwa aln -t 2 Ref paired_data_(varant name)1.fq > (variant name)_sequence1.sai
```
```bash
bwa aln -t 2 Ref paired_data_(varant name)2.fq > (variant name)_sequence2.sai
```
Generate SAM files for the aligned reads: 
```bash
bwa sampe Ref (variant name)_sequence1.sai (variant name)_sequence2.sai (variant name)_sequence1.fq. (variant name)_sequence2.fq > (variant name)_ sequence12_pe.sam
```
Convert the SAM file to BAM file using [SamTools]( http://www.htslib.org/doc/). 
```bash
sudo apt-get install bwa
```
```bash
samtools view -b (variant name)_sequence12_pe.sam > (variant name).bam
```
Rebuild with SamTools to sort genome: 
```bash
samtools sort -o (variant name)_sorted.bam (variant name).bam
```
Use freebayes to generate VCF files vs. reference for each variant: 
```bash
sudo apt-get install freebayes
```
```bash
freebayes -f GCF_009858895.2_ASM985889v3_genomic.fna. (variant name)_sorted.bam > (variant name).vcf
```

Transfer the VCF files to R for statistical analysis. From the VCF files, extract the allele and locus identity were extract for each genetic variation from the 5 source genomes (4 variants and the reference) in a excel file. The phenotypes (transmissibility) of the variants were collected from primary literature. Example excel table for the first variant locus:  

Source | Allele | Locus | Phenotype  
--- | --- | --- | --- 
Ref | G | 174| 1.270 
UK  | G| 174 | 2.108 
Brazil | G | 174 | 3.874 
California | G | 174 | 1.541
SA| T | 174 | 1.980 



