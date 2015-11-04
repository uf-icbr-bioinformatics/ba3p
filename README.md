# ba3p
### Bacterial alignment, assembly, and annotation pipelines

**Ba3p** is a collection of pipeline scripts to perform high-throughput analysis of NGS data from bacterial genomes. It currently includes the following tools:

- [snpcall](#snpcall)

In addition, this package provides a few utilities that are used by ba3p tools, but can also be used independently if desired:

- VCFmerger.py
- sam-remap.py
- vcf-remap.py

## Requirements
Ba3p is designed to run in an HPC environment, and requires the **submit** and **actor** packages.  Please contact the author for details.

## Basic usage

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
