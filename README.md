# ba3p - Bacterial alignment, assembly, and annotation

**Ba3p** is a collection of pipelines and scripts to perform high-throughput analysis of NGS data from bacterial genomes. It currently includes the following tools:

- [snpcall](#user-content-snpcall)
- [denovo](#user-content-denovo)

In addition, this package provides a few utilities that are used by ba3p tools, but can also be used independently if desired:

- [swipt.py](#user-content-swiptpy)
- [VCFmerger.py](#user-content-vcfmergerpy)
- [sam-remap.py](#user-content-samremappy)
- [vcf-remap.py](#user-content-vcfremappy)
- [writeModelXML.py](#user-content-writemodelxmlpy)

### Requirements
Ba3p is designed to run in an HPC environment, and requires the **submit** and **actor** packages.  Please contact the author for details.

### Basic usage

All ba3p tools are invoked as follows:

```
$ ba3p <toolName> <configurationFile>
```

You can also enter:

```
$ ba3p help <toolName>
```

to access the interactive manual for each tool.

The rest of this document describes the currently available tools and explains their usage.

## Snpcall

The **snpcall** tool performs the following processing steps on one or more pairs of (paired-end) FASTQ files:

- Trim short reads and align them to the reference using bowtie2;
- Remove duplicates with GATK;
- SNP calling with freebayes;
- Optional: SNP annotation using snpEff.

The above steps can be performed in parallel for any number of samples. At the end,
you have the option to combine the VCF files produced by freebayes for all the
samples, generating a set of multi-fasta files suitable for phylogenetic analysis.

### Setup
The snpcall tool assumes you have the following files: 

- Two (paired-end) fastq files for each run. The fastq files may be compressed 
  with gzip. [Note: single-end support is coming soon]
- The reference sequence of the genome being sequenced.
- The bowtie2 index for the reference sequence.

Optionally: 

- The name of the snpEff database associated with the genome.
- A mapping file, converting chromosome names used in the reference sequence
  to those used by snpEff, if they are different.

### Configuration
Snpcall requires a configuration file listing the parameters
of the desired analysis. The configuration file is a plain-text file divided
into sections. Each section starts with a label surrounded by [ ], and contains
assignments of values to variables in the form: variable = value.

The first section must be called **[General]**, and must contain the following
variables:

```conf
title = (name of your analysis)
reference = (path to reference sequence in fasta format)
genome = (path to bowtie2 index of genome)
samples = (names of samples separated by commas).
```

In addition, it may contain:

```conf
snpeffdb = (name of snpEff database)
mapfile = (mapping between reference and snpEff chromosome names).
```
The remainder of the configuration file is composed of one section
for each sample, labeled with its name. For example, if you defined:

```conf
samples = smpl1,smpl2
```

the corresponding sections should be labeled **[smpl1]**, **[smpl2]**.

Each section should indicate the name of the two paired-end fastq files:
```conf
fastq1 = (file with left-side reads)
fastq2 = (file with righ-side reads)
```
Please note that all pathnames should be either absolute (ie, starting
with a /) or relative to the directory you are calling the script from.

### Output
When snpcall starts it will create a directory with the 
name indicated by the 'title' variable, and will write all its output 
files to that directory.

At the end, the tool will create a zip file with the contents of the
output directory, excluding temporary files and other unnecessary files.
The zip file will also include an HTML file (called **index.html**) containing
a full report of the analysis that was performed, and links to all output
files.

The easiest way to access the pipeline results is to download the zip file 
to your computer and uncompress it, then open the index.html file with 
any web browser.

### Example
This is an example configuration file for two pairs of fastq files
from the human genome. You can copy and paste the following text to a new file
and use it as a template to create your own.

```conf
[General]
title = testrun
reference = /lfs/scratch/bio/reference/samtools/hg19.fa
genome = /lfs/scratch/bio/reference/bowtie2/hg19
snpeffdb = GRCh38.76
samples = sample1,sample2

[sample1]
fastq1=TESTRUN_S1_L001_R1.fastq.gz
fastq1=TESTRUN_S1_L001_R2.fastq.gz

[sample2]
fastq1=TESTRUN_S2_L001_R1.fastq.gz
fastq1=TESTRUN_S2_L001_R2.fastq.gz
```

## Denovo

The **denovo** tool performs the following processing steps on one or more pairs of (paired-end) FASTQ files:

- Quality control on raw reads using **FastQC**;
- Trimming of raw reads using **sickle**;
- Quality control on trimmed reads using **FastQC**;
- Optimized de-novo assembly of trimmed reads using **VelvetOptimizer**.

### Setup
The denovo tool assumes you have the following files:

- Two (paired-end) fastq files for each run. The fastq files may be compressed 
  with gzip. [Note: single-end support is coming soon]

### Configuration
Denovo requires a configuration file listing the parameters
of the desired analysis. The configuration file is a plain-text file divided
into sections. Each section starts with a label surrounded by [ ], and contains
assignments of values to variables in the form: variable = value.

The first section must be called **[General]**, and must contain the following
variables:

```conf
title = (name of your analysis)
samples = (names of samples separated by commas).
```
The remainder of the configuration file is composed of one section
for each sample, labeled with its name. For example, if you defined:

```conf
samples = smpl1,smpl2
```

the corresponding sections should be labeled **[smpl1]**, **[smpl2]**.

Each section should indicate the name of the two paired-end fastq files:
```conf
fastq1 = (file with left-side reads)
fastq2 = (file with righ-side reads)
```
Please note that all pathnames should be either absolute (ie, starting
with a /) or relative to the directory you are calling the script from.

### Output
When denovo starts it will create a directory with the 
name indicated by the 'title' variable, and will write all its output 
files to that directory.

At the end, the tool will create a zip file with the contents of the
output directory, excluding temporary files and other unnecessary files.
The zip file will also include an HTML file (called **index.html**) containing
a full report of the analysis that was performed, and links to all output
files.

The easiest way to access the pipeline results is to download the zip file 
to your computer and uncompress it, then open the index.html file with 
any web browser.

## swipt.py
**swipt.py** performs sliding-window analysis of recombination using the Phi test, as implemented by T. Bruen and D. Bryant (see the **Phi test** section in [this page](http://www.maths.otago.ac.nz/~dbryant/software.html)). Given a multiple sequence alignment in FASTA format as input, swipt applies the Phi test to a sliding subregion of the alignment, and reports the resulting P-value for each region to a tab-delimited file (or to standard output).

### Syntax
swipt.py is invoked as follows:
```
swipt.py [options] alignment.fa [outfile]
```

The following command-line options are available:

Option  | Description
--------|------------
-p PATH | Set the path to the Phi executable
-w W    | Set window size to W (default: 200)
-s S    | Set window step to S (default: 40)
-t T    | Set P-value threshold to T. If this option is specified, only regions with a Phi P-value below this limit will be reported in the output.
-phip N | Set Phi permutations to N (Passed to Phi as the -p option)
-phio   | Output extended Phi results (Passed to Phi as the -o option)

### Output
The output is a tab-delimited file containing either three or five columns. The first two columns contain the start and end position of each region, while the third column contains the P-value from the Phi test. If *-phio* is specified (corresponding to the -o option in the Phi program), the fourth column will contain the P-value from the *NSS* test, and the fifth column will contain the P-value from the *Max Chi<sup>2</sup>* test. Please note that the P-value threshold specified with the -t option only applies to the Phi test result (third column).

## VCFmerger.py
**VCFmerger.py** is a program to combine multiple VCF files (each one possibly containing data for one or more samples) into pseudo-FASTA files. A pseudo-fasta file is a FASTA file containing multiple sequences, one for each sample, and each sequence is composed of the nucleotides found in variant positions only. For example, if the input to VCFmerger consists of three VCF files (containing a single sample each) with 10 variant positions in total, the output file will be a pseudo-FASTA file containing three sequences of 10nt each.

### Syntax
VCFmerger.py is invoked as follows:
```
VCFmerger.py [options] vcffiles ...
```
The following command-line options are available (each option has both a short and a long form):

Short | Long | Description
------|------|------------
  -m  |--missingChar C  | Use character C to represent unknown alleles (otherwise the reference allele will be used).
  -o  |--outPattern P   | Use string P as the name of the output file(s). If P contains the sequence {}, it will be replaced by each chromosome name, and a separate output file will be created for each chromosome. Otherwise, output will be written to a single file.
  -b  |--bothAlleles    | If specified, the program will print two bases for each SNP, according to the genotype (AA, AB, or BB).
  -r  |--removeMissing  | If specified, the program will only print SNPs that appear in all the VCF files.
  -q  |--quality Q      | Use Q as the quality threshold. If not specified, the quality threshold will be read from the VCF file, if present, otherwise it will be set to the default value of 30.
  -R  |--report R       | Write the list of output files produced to a report file R in tab-delimited format. Each line contains three columns: chromosome, filename, number of SNPs.
  -i  |--impact I       | If snpEff annotations are present, only retain SNPs with a putative impact I or higher. I should be one of HIGH, MODERATE, LOW, or MODIFIER.
  -a|--annotation A[,A...]  | If snpEff annotations are present, only retain SNPs having an annotation matching one or more of the supplied annotations A. See http://snpeff.sourceforge.net/VCFannotationformat_v1.0.pdf for a list of all valid annotations, or use the -l option.
  -s|--saveSNPs S           | Save SNPs to tab-delimited file S.
  -l|--listEff              | List all valid values for snpEff impact and annotation. When the long form is used, also prints a short description of each annotation (if available).

## sam-remap.py
**sam-remap.py** converts chromosome names in a SAM file. It is invoked as follows:
```
sam-remap.py mapfile
```

The file *mapfile* should be tab-delimited file with two columns: the first column
contains chromosome names as they currently appear in the SAM file, and the
second column contains the new chromosome names. For example, to rename chromosomes
from roman numerals to standard form:

```
chrI	chr1
chrII	chr2
...
```

This command is meant to work as a filter: it reads a SAM file from standard
input and writes a new SAM file to standard output. For example:

```
> cat existing.sam | sam-remap.py map.txt > converted.sam
```

To convert chromosome names in a BAM file use something like the following:

```
> samtools view -h existing.bam | sam-remap.py map.txt | samtools view -b -o converted.bam -
```

## vcf-remap.py
**vcf-remap.py** converts chromosome names in a VCF file. It is invoked as follows:
```
vcf-remap.py mapfile
```

The file *mapfile* should be tab-delimited file with two columns: the first column
contains chromosome names as they currently appear in the VCF file, and the
second column contains the new chromosome names. For example, to rename chromosomes
from roman numerals to standard form:

```
chrI	chr1
chrII	chr2
...
```

This command is meant to work as a filter: it reads a VCF file from standard
input and writes a new VCF file to standard output. For example:

```
  > cat existing.vcf | vcf-remap.py map.txt > converted.vcf
```

## writeModelXML.py
**writeModelXML.py** is a tool to generate the linear model section of a BEAST XML file.
It writes an XML fragment containing the three sections **generalDataType**, **glmSubstitutionModel**, and 
**markovJumpsTreeLikelihood**. This fragment should then be inserted into the overall BEAST file.

The program can be used in two modes, *standard* or *generate*.

### Standard mode

Usage:
```
writeModelXML.py [-i] infile [outfile]
```
The input file `infile' should be in tab-delimited format and should contain location names, 
one per line. Additional columns represent descriptors. The first line of the input file 
should contain descriptor names. For example:
```
Location   Desc1   Desc2   Desc3
LOC1       0.13    1.32    -0.3
LOC2       -1.14   0.21    -2.0
LOC3       2.11    0.67    -0.55
```
XML output will be written to `outfile' if specified, or to standard output.

The -i option causes an inverse descriptor to be added for each descriptors specified
in the input file.

### Generate mode

Usage: 
```
writeModelXML.py -g outfile nlocs ndescs
```

If invoked with the -g option, the program will create an empty locations file `outfile' for 
*nlocs* locations and *ndescs* descriptors. This file can then be used as an input file for
standard mode, after filling in the descriptor columns.
