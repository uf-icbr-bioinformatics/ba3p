# (c) 2015, A. Riva, DiBiG

import csv
import os.path

def message(string, *args):
    # Write `string' to standard error. `args' are 
    # inserted into `string' with format.
    sys.stderr.write(string.format(*args))
    sys.stderr.write("\n")

def linkify(url, name=False):
    if not name:
        name = url
    return "<A href='{}'>{}</A>".format(url, name)

def fixPath(path):
    """Add ../ in front of `path' unless it is absolute."""
    if path[0] == "/":
        return path
    else:
        return "../" + path

class ba3pMSC():
    valid = True                # Set to False in case of configuration errors
    title = "snpcall"           # Title of run
    runs = []                   # List of runs
    nruns = 0                   # Number of runs
    reference = None            # Genome reference sequence (FASTA)
    genome = None               # Bowtie2 index of genome
    mapfile = None              # For chromosome renaming
    snpeffdb = None             # For snpEff annotations

    def __init__(self, ACT):
        global message
        conffile = ACT.Arguments[0]
        conf = ACT.loadConfiguration(conffile)
        self.title = ACT.getConf("title")
        self.reference = ACT.getConf("reference")
        self.genome = ACT.getConf("genome")
        self.mapfile = ACT.getConf("mapfile")
        self.snpeffdb = ACT.getConf("snpeffdb")
        self.singleVCF = ACT.getConf("singleVCF")

        # For singleVCF mode
        self.vcfsnps = None
        self.seVCF = None
        self.snpeffVCF = None
        self.csvreport = None
        self.htmlreport = None

        if ACT.missingOrStale(self.genome, warn=True):
            valid = False
        if ACT.missingOrStale(self.reference + ".1.bt2", warn=True):
            valid = False
        if self.mapfile and ACT.missingOrStale(self.mapfile, warn=True):
            valid = False

        runnames = ACT.getConf("samples").split(",")
        runnames = [ s.strip() for s in runnames ]

        for name in runnames:
            run = {'name': name,
                   'fastq1': ACT.getConf("fastq1", name),
                   'fastq2': ACT.getConf("fastq2", name)}

            if run['fastq1'] == None or ACT.missingOrStale(run['fastq1']):
                message("Configuration error: `fastq1' undefined for sample `{}'".format(name))
                valid = False
            if run['fastq2'] == None or ACT.missingOrStale(run['fastq2']):
                message("Configuration error: `fastq2' undefined for sample `{}'".format(name))
                valid = False

            self.runs.append(run)
            self.nruns += 1

def dumpmsc(msc):
    global message
    message("Title: {}", msc.title)
    message("Reference: {}", msc.reference)
    message("Genome index: {}", msc.genome)
    message("Mapfile: {}", msc.mapfile)
    message("SnpEffDB: {}", msc.snpeffdb)
    message("{} samples:", msc.nruns)
    for run in msc.runs:
        message("  Name: {}", run['name'])
        message("  Fastq1: {}", run['fastq1'])
        message("  Fastq2: {}", run['fastq2'])
        message("")

MSC = ba3pMSC(ACT)

if not MSC.valid:
    message("Script cannot be executed due to configuration errors.")
    exit
dumpmsc(MSC)

# Script definition 

ACT.script(MSC.title, "BA3P - SNP calling and annotation", "BA3P")
ACT.begin(timestamp=False)

ACT.scene(1, "General configuration")
ACT.reportf("""Sample name: <b>{}</b><br>
Reference genome: <b>{}</b><br>
Bowtie2 index: <B>{}</b><br>
""".format(MSC.title, MSC.reference, MSC.genome))
if MSC.mapfile != None:
    ACT.reportf("Chromosome mapping: <b>{}</b><br>".format(MSC.mapfile))
if MSC.snpeffdb != None:
    ACT.reportf("snpEff database: <b>{}</b><br>".format(MSC.snpeffdb))
ACT.reportf("Samples and input files:<br>")
ACT.table([ [r['name'], r['fastq1'], r['fastq2'] ] for r in MSC.runs],
          header=["Sample", "Left reads", "Right reads"],
          align="HLL")

MSC.reference = fixPath(MSC.reference)
MSC.genome = fixPath(MSC.genome)
if MSC.mapfile != None:
    MSC.mapfile = fixPath(MSC.mapfile)

#
# Before we start, check if reference file is indexed
#

fai = MSC.reference + ".fai" 
refdict = ACT.setFileExt(MSC.reference, ".dict")

if ACT.missingOrStale(fai, MSC.reference):
    print "Reference index {} does not exist or is out of date, creating it.".format(fai)
    ACT.shell("module load samtools; samtools faidx {}", MSC.reference)

if ACT.missingOrStale(refdict, MSC.reference):
    print "Reference dictionary {} does not exist, creating it.".format(refdict)
    ACT.submit("picard.qsub CreateSequenceDictionary R={} O={}".format(MSC.reference, refdict), done="dict.done")
    ACT.wait("dict.done")

# Ensure we don't have old .done files lying around
ACT.shell("rm -f *.done")

#
# First of all run sickle and bowtie on fastq files
#

for run in MSC.runs:
    name = run['name']
    fastq1 = fixPath(run['fastq1'])
    fastq2 = fixPath(run['fastq2'])

    in1 = os.path.basename(fastq1)
    in2 = os.path.basename(fastq2)
    run['sickle1'] = ACT.setFileExt(in1, ".sickle.fastq", remove=[".fastq", ".gz"])
    run['sickle2'] = ACT.setFileExt(in2, ".sickle.fastq", remove=[".fastq", ".gz"])
    run['bamfile'] = run['name'] + ".bam"
    this = ACT.submit("sickle.qsub {} {} {} {}".format(fastq1, fastq2, run['sickle1'], run['sickle2']))
     # print "Submitted: " + this
    ACT.submit("bowtie2-pe.qsub {} {} {} {}".format(MSC.genome, run['sickle1'], run['sickle2'], run['bamfile']), done="bowtie2.@.done", after=this)
ACT.wait(('bowtie2.@.done', MSC.nruns))

#
# Once all the bowties are done, we can collect stats (number of aligned reads),
# and then proceed with the GATK pipeline, followed by SNP calling.
#

for run in MSC.runs:
    ACT.submit("gzip.qsub {}".format(run['sickle1']))
    ACT.submit("gzip.qsub {}".format(run['sickle2']))

    name = run['name']
    bamfile = run['bamfile']
    run['alignedReads'] = ACT.shell("module load bamtools; bamtools count -in {}".format(bamfile))
    
    run['picard'] = ACT.setFileExt(bamfile, ".picard.bam")
    job = ACT.submit("picard.qsub AddOrReplaceReadGroups I={} O={} RGID=ID_{} RGLB=LB_{} RGPL=ILLUMINA RGPU=PU_{} RGSM=SM_{}".format(bamfile, run['picard'], name, name, name, name))

    run['dupmark'] = ACT.setFileExt(bamfile, ".dupmark.bam")
    run['metrics'] = name + ".metrics.txt"
    job = ACT.submit("picard.qsub MarkDuplicates I={} O={} M={} REMOVE_DUPLICATES=true ASSUME_SORTED=true".format(run['picard'], run['dupmark'], run['metrics']), 
                     after=job)

    run['realigned'] = name + ".realigned.bam"
    job = ACT.submit("gatk-realigner.qsub {} {} {}".format(MSC.reference, run['dupmark'], run['realigned']), after=job)

# Fix mates

    run['fixmates'] = name + ".fixmates.bam"
    ACT.submit("picard.qsub FixMateInformation I={} O={} SORT_ORDER=coordinate".format(run['realigned'], run['fixmates']), after=job, done="picard.@.done")
ACT.wait(("picard.@.done", MSC.nruns))

# And now the main star...

if MSC.singleVCF:
    # Start a single freebayes run for all samples
    bamfiles = ",".join([ run['fixmates'] for run in MSC.runs ])
    MSC.singleVCF = "all.vcf"
    ACT.submit("freebayes.qsub {} {} {}".format(MSC.reference, bamfiles, MSC.singleVCF), done="freebayes.done")
    ACT.wait("freebayes.done")
else:
    # Start a separate freebayes run for each sample
    for run in MSC.runs:
        run['vcf'] = name + ".vcf"
        ACT.submit("freebayes.qsub {} {} {}".format(MSC.reference, run['fixmates'], run['vcf']), done="freebayes.@.done")
    ACT.wait(("freebayes.@.done", MSC.nruns))

# Let's collect some stats
for run in MSC.runs:
    run['fixedReads'] = ACT.shell("module load bamtools; bamtools count -in {}".format(run['fixmates']))
    if not MSC.singleVCF:
        run['vcfsnps'] = ACT.shell("grep -c -v ^\# {}".format(run['vcf']))

if MSC.singleVCF:
    MSC.vcfsnps = ACT.shell("grep -c -v ^\# {}".format(MSC.singleVCF))

ACT.scene(2, "Alignment")
ACT.reportf("""Reads were aligned to the reference genome using <b>bowtie2</b>, after trimming with <b>sickle</b>.
Aligned reads were processed with the <b>GATK</b> pipeline, to remove duplicates and fix mates.
The following table shows the original number of aligned reads and the number of reads after GATK processing for each sample.
Each sample name is linked to the corresponding BAM file.""")
ACT.table([ ["<a href='{}'>{}</a>".format(run['fixmates'], run['name']), run['alignedReads'], run['fixedReads'] ] for run in MSC.runs], 
          header=['Sample', 'Aligned reads', 'Post-GATK reads'], align="HRR")

ACT.scene(3, "SNP calling")
if MSC.singleVCF:
    ACT.reportf("""SNPs were called using <b>freebayes</b> in single-VCF mode (all BAM files were combined
together for this step. A total of {} SNPs were identified. The VCF file can be downloaded using the link below.""", MSC.vcfsnps)
    ACT.file(MSC.singleVCF, description="VCF file containing all identified SNPs.")
else:
    ACT.reportf("""SNPs were called using <b>freebayes</b>. The following table shows the number of SNPs
identified in each sample. Each sample name is linked to the corresponding VCF file.""")
    ACT.table([ ["<a href='{}'>{}</a>".format(r['vcf'], r['name']), r['vcfsnps'] ] for r in MSC.runs],
              header=['Sample', 'SNPs'], align='HR')

#
# If requested, call snpEff
#

if MSC.snpeffdb != None:
    if MSC.singleVCF:
        snpeffIn = MSC.singleVCF
        # Rename chromosomes if requested
        if MSC.mapfile != None:
            MSC.seVCF = "all.se.vcf"
            snpeffIn = MSC.seVCF
            ACT.shell("module load dibig_ba3p; vcf-remap.py {} < {} > {}".format(MSC.mapfile,  MSC.singleVCF, MSC.seVCF))
        ACT.mkdir("allsnps")
        MSC.snpeffVCF = "all.snpeff.vcf"
        MSC.csvreport = "allsnps/allsnps.csv"
        MSC.htmlreport = "allsnps/allsnps.html"
        ACT.submit("snpeff.qsub {} {} {} csv={} html={}".format(snpeffIn, MSC.snpeffVCF, MSC.snpeffdb, MSC.csvreport, MSC.htmlreport), done="snpeff.done")
        ACT.wait("snpeff.done")
    else:

        # One VCF run per sample
        for run in MSC.runs:
            name = run['name']
            vcf = run['vcf']
            snpeffIn = vcf
            # Rename chromosomes if requested
            if MSC.mapfile != None:
                run['vcfse'] = vcfse = name + ".se.vcf"
                snpeffIn = vcfse
                ACT.shell("module load dibig_ba3p; vcf-remap.py {} < {} > {}".format(MSC.mapfile,  vcf, vcfse))
            
            ACT.mkdir(name)
            run['snpeff'] = name + ".snpeff.vcf"
            run['csvreport'] = name + "/" + name + ".csv"
            run['htmlreport'] = name + "/" + name + ".html"
            ACT.submit("snpeff.qsub {} {} {} csv={} html={}".format(snpeffIn, run['snpeff'], MSC.snpeffdb, run['csvreport'], run['htmlreport']), done="snpeff.@.done")
        ACT.wait(("snpeff.@.done", MSC.nruns))

lastscene = 4
if MSC.snpeffdb != None:
    lastscene = 5
    ACT.scene(4, "SNP functional annotation")
    if MSC.singleVCF:
        ACT.reportf("""SNPs were annotated <b>snpEff</b> using the <b>{}</b> database. The following table provides links to the annotated VCF file
and to the snpEff reports (in CSV and HTML formats).""".format(MSC.snpeffdb))
        ACT.table([[linkify(MSC.snpeffVCF), linkify(MSC.csvreport), linkify(MSC.htmlreport)]],
                  header=['VCF', 'CSV report', 'HTML report'],
                  align="LLL")
    else:
        ACT.reportf("""SNPs were annotated through <b>snpEff</b> using the <b>{}</b> database. The following table provides links to the annotated VCF file
and to the snpEff reports (in CSV and HTML formats) for each sample.""".format(MSC.snpeffdb))
        ACT.table([ [r['name'], linkify(r['snpeff']), linkify(r['csvreport']), linkify(r['htmlreport']) ] for r in MSC.runs],
                  header=['Sample', 'VCF', 'CSV report', 'HTML report'],
                  align="HLLL")

#
# Finally combine all VCFs
#

allvcfs = ""
report = "merged-snps.csv"

if MSC.singleVCF:
    if MSC.snpeffdb == None:
        invcf = MSC.singleVCF
    else:
        invcf = MSC.snpeffVCF
    ACT.shell("module load dibig_ba3p; VCFmerger.py -R {} {}".format(report, invcf))
else:
    for run in MSC.runs:
        if 'vcfse' in run:
            allvcfs = allvcfs + run['vcfse'] + " "
        else:
            allvcfs = allvcfs + run['vcf'] + " "
    ACT.shell("module load dibig_ba3p; VCFmerger.py -R {} {}".format(report, allvcfs))

ACT.scene(lastscene, "Combined SNP files")
ACT.reportf("""The VCF files for all samples were combined into a set of multi-fasta files (one for each chromosome).
The following table provides a link to the FASTA file for each chromosome, with the number of SNPs it contains.""")

# This should be turned into a function...
data = []
first = True
f = open(report, "r")
rows = csv.reader(f, delimiter='\t')
for row in rows:
    if first:
        first = False
    else:
        row[1] = linkify(row[1])
        data.append(row)
ACT.table(data, header=['Chromosome', 'FASTA', 'SNPs'], align="HLR")
f.close()

# Should do some cleanup at the end...

toDelete = ""
for run in MSC.runs:
    toDelete += " " + run['picard']
    toDelete += " " + run['dupmark']
    toDelete += " " + run['realigned']
ACT.execute("rm" + toDelete)

