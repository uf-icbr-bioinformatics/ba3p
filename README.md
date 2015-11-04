# ba3p
### Bacterial alignment, assembly, and annotation pipelines

**Ba3p** is a collection of pipeline scripts to perform high-throughput analysis of NGS data from bacterial genomes. It currently includes the following tools:

- snpcall

## Requirements
Ba3p is designed to run in an HPC environment, and requires the **submit** and **actor** packages.  Please contact the author for details.

## Basic usage

All ba3p tools are invoked as follows:

```
$ ba3p <toolName> <configurationFile>
```

You can also enter:

```bash
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
The snpcall script assumes you have the following files: 

- Two (paired-end) fastq files for each run. The fastq files may be compressed 
  with gzip. [Note: single-end support is coming soon]
- The reference sequence of the genome being sequenced.
- The bowtie2 index for the reference sequence.

Optionally: 

- The name of the snpEff database associated with the genome.
- A mapping file, converting chromosome names used in the reference sequence
  to those used by snpEff, if they are different.
