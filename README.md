# ba3p - Bacterial alignment, assembly, and annotation

**Ba3p** is a collection of pipelines and scripts to perform high-throughput analysis of NGS data from bacterial genomes. It currently includes the following tools:

- [snpcall](#user-content-snpcall)
- denovo

In addition, this package provides a few utilities that are used by ba3p tools, but can also be used independently if desired:

- swipt.py
- VCFmerger.py
- sam-remap.py
- vcf-remap.py

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
  -o  |--outPattern P   | Use string P as the name of the output file(s). If P contains the sequence {}, it will be replaced by each chromosome name.
  -b  |--bothAlleles    | If specified, the program will print two bases for each SNP, according to the genotype (AA, AB, or BB).
  -r  |--removeMissing  | If specified, the program will only print SNPs that appear in all the VCF files.
  -q  |--quality Q      | Use Q as the quality threshold. If not specified, the quality threshold will be read from the VCF file, if present, otherwise it will be set to the default value of 30.
  -R  |--report R       | Write the list of output files produced to a report file R in tab-delimited format. Each line contains three columns: chromosome, filename, number of SNPs.
  -i  |--impact I       | If snpEff annotations are present, only retain SNPs with a putative impact I or higher. I should be one of HIGH, MODERATE, LOW, or MODIFIER.
  -a|--annotation A[,A...]  | If snpEff annotations are present, only retain SNPs having an annotation matching one or more of the supplied annotations A. See http://snpeff.sourceforge.net/VCFannotationformat_v1.0.pdf for a list of all valid annotations, or use the -l option.
  -s|--saveSNPs S           | Save SNPs to tab-delimited file S.
  -l|--listEff              | List all valid values for snpEff impact and annotation. When the long form is used, also prints a short description of each annotation (if available).
