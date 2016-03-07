# MultiSampleActor

# (c) 2015, A. Riva, DiBiG
# Multi-sample SNP calling pipeline, v2.0

# ToDo: use trimmomatic instead of sickle

import sys
import os.path

from SampleCollection import SampleCollection

def linkify(url, name=False):
    if not name:
        name = url
    return "<A href='{}'>{}</A>".format(url, name)

# Methods we add to the ACT object

def dump(self):
    self.message("Title: {}", self.title)
    self.message("Reference: {}", self.reference)
    self.message("Genome index: {}", self.btidx)
    self.message("Mapfile: {}", self.mapfile)
    self.message("SnpEffDB: {}", self.snpeffdb)
    self.sc.showSamples()

def sickleBowtie(self):
    for rs in self.sc.readsets:
        fastq1 = ACT.fixPath(rs['left'])
        fastq2 = ACT.fixPath(rs['right'])

        in1 = os.path.basename(fastq1)
        in2 = os.path.basename(fastq2)
        rs['sickle1'] = self.setFileExt(in1, ".sickle.fastq.gz", remove=[".fastq", ".gz"])
        rs['sickle2'] = self.setFileExt(in2, ".sickle.fastq.gz", remove=[".fastq", ".gz"])
        rs["fqcdir1"]   = in1 + ".fqc"
        rs["fqcdir2"]   = in2 + ".fqc"
        self.mkdir(rs["fqcdir1"])
        self.mkdir(rs["fqcdir2"])
        rs['bamfile'] = self.exclude(rs['name'] + ".bam")
        if self.missingOrStale(rs['sickle1']) or self.missingOrStale(rs['sickle2']): # Don't redo sickle if already done
            this = self.submit("sickle.qsub pe {} {} {} {}".format(fastq1, fastq2, rs['sickle1'], rs['sickle2']))
            self.submit("fastqc.qsub {} {}".format(rs["sickle1"], rs["fqcdir1"]), after=this, done="fqc2.@.done")
            self.submit("fastqc.qsub {} {}".format(rs["sickle2"], rs["fqcdir2"]), after=this, done="fqc2.@.done")
            self.submit("bowtie2-pe.qsub {} {} {} {}".format(self.btidx, rs['sickle1'], rs['sickle2'], rs['bamfile']), done="bowtie2.@.done", after=this)
        else:
            self.submit("bowtie2-pe.qsub {} {} {} {}".format(self.btidx, rs['sickle1'], rs['sickle2'], rs['bamfile']), done="bowtie2.@.done")
    self.wait(('bowtie2.@.done', self.sc.nreadsets))

def postBowtie(self):
    """Check that the bowtie output files exist and count the number of aligned reads in each.
NOTE: this should be replaced with countBAMs from rnaseq.py."""
    for rs in self.sc.readsets:
        bamfile = rs['bamfile']
        if ACT.checkFile(bamfile, step="bowtie2"):
            rs['alignedReads'] = self.shell("module load bamtools; bamtools count -in {}".format(bamfile))

def doVCFsingleOld(self):
    # Start a single freebayes run for all samples
    bamfiles = ",".join([ smp['fixmates'] for smp in self.sc.readsets ])
    self.singleVCF = "all.vcf"
    self.submit("freebayes.qsub {} {} {}".format(self.reference, bamfiles, self.singleVCF), done="freebayes.done")
    self.wait("freebayes.done")
    self.vcsnps = self.shell("grep -c -v ^\# {}".format(self.singleVCF))

def doVCFsingle(self):
    self.singleVCF = "all.vcf"
    with open("VCFLIST", "w") as v:
        for rs in self.sc.readsets:
            rs['vcf'] = self.exclude(rs['name'] + ".g.vcf")
            v.write(rs['vcf'] + "\n")
            self.submit("gatk-snpcall-gvcf.qsub {} {} {}".format(self.reference, rs['bamfile'], rs['vcf']), done="freebayes.@.done")
    self.wait(("freebayes.@.done", self.sc.nreadsets))
    self.submit("gvcf.qsub {} {} {}".format(self.reference, "VCFLIST", self.singleVCF), done="gvcf.done")
    self.wait("gvcf.done")
    for smp in self.sc.readsets:
        smp['vcfsnps'] = 0
        self.submit("gzip.qsub {}".format(smp['vcf']))
    self.vcsnps = self.shell("grep -c -v ^\# {}".format(self.singleVCF))

def doVCFmulti(self):
    for smp in self.sc.readsets:
        smp['vcf'] = smp['name'] + ".vcf"
        self.submit("freebayes.qsub {} {} {}".format(self.reference, smp['fixmates'], smp['vcf']), done="freebayes.@.done")
    self.wait(("freebayes.@.done", self.sc.nreadsets))
    for smp in self.sc.readsets:
        smp['vcfsnps'] = self.shell("grep -c -v ^\# {}".format(smp['vcf']))

def snpeffSingle(self):
    snpeffIn = self.singleVCF
    # Rename chromosomes if requested
    if self.mapfile != None:
        self.seVCF = "all.se.vcf"
        snpeffIn = self.seVCF
        self.shell("module load dibig_ba3p; vcf-remap.py {} < {} > {}".format(self.mapfile,  self.singleVCF, self.seVCF))
    self.mkdir("allsnps")
    self.snpeffVCF = "all.snpeff.vcf"
    self.csvreport = "allsnps/allsnps.csv"
    self.htmlreport = "allsnps/allsnps.html"
    self.submit("snpeff.qsub {} {} {} csv={} html={} ver={} dir=P".format(snpeffIn, self.snpeffVCF, self.snpeffdb, self.csvreport, self.htmlreport, self.snpeffver), done="snpeff.done")
    self.wait("snpeff.done")

def snpeffMulti(self):
    for smp in self.sc.readsets:
        name = smp['name']
        vcf = smp['vcf']
        snpeffIn = vcf
        # Rename chromosomes if requested
        if self.mapfile != None:
            smp['vcfse'] = vcfse = name + ".se.vcf"
            snpeffIn = vcfse
            self.shell("module load dibig_ba3p; vcf-remap.py {} < {} > {}".format(self.mapfile,  vcf, vcfse))
            
        self.mkdir(name)
        smp['snpeff'] = name + ".snpeff.vcf"
        smp['csvreport'] = name + "/" + name + ".csv"
        smp['htmlreport'] = name + "/" + name + ".html"
        self.submit("snpeff.qsub {} {} {} csv={} html={} ver={} dir=P".format(snpeffIn, smp['snpeff'], self.snpeffdb, smp['csvreport'], smp['htmlreport'], self.snpeffver), done="snpeff.@.done")
    self.wait(("snpeff.@.done", self.sc.nreadsets))

def doVCFmergeSingle(self):
    if self.snpeffdb == None:
        invcf = self.singleVCF
    else:
        invcf = self.snpeffVCF
    self.submit('generic.qsub VCFmerger.py -R {} {} module:dibig_ba3p'.format(self.report, invcf), done="merger.done")
    self.wait("merger.done")

def doVCFmergeMulti(self):
    allvcfs = ""
    for smp in self.sc.readsets:
        if 'vcfse' in smp:
            allvcfs = allvcfs + smp['vcfse'] + " "
        else:
            allvcfs = allvcfs + smp['vcf'] + " "
    self.submit('generic.qsub VCFmerger.py -R {} {} module:dibig_ba3p'.format(self.report, allvcfs), done="merger.done")
    self.wait("merger.done")

def readReport(self):
    global linkify
    data = []
    first = True
    with open(self.report, "r") as f:
        for row in f:
            if first:
                first = False
            else:
                parsed = row.split('\t')
                parsed[1] = linkify(parsed[1])
                data.append(parsed)
    return data

### Script code starts here

## Initialization
ACT.loadConfiguration(ACT.Arguments[0])
SC = SampleCollection(ACT.Conf)
ACT.sc = SC

VERIFYFILES = True

## Global configuration
ACT.title = ACT.getConf("title")
ACT.reference = ACT.checkPath(ACT.getConf("reference"))
ACT.checkPath(ACT.getConf("btidx") + ".1.bt2")
ACT.btidx = ACT.fixPath(ACT.getConf("btidx"))
ACT.mapfile = ACT.checkPath(ACT.getConf("mapfile"))
ACT.snpeffdb = ACT.getConf("snpeffdb")
ACT.snpeffver = ACT.getConf("snpeffver", default="4.1")
ACT.singleVCF = ACT.getConf("singleVCF")

## For singleVCF mode
ACT.vcfsnps = None
ACT.seVCF = None
ACT.snpeffVCF = None
ACT.csvreport = None
ACT.htmlreport = None
ACT.report = "merged-snps.csv"

## Check that all fastq files exist
for r in SC.readsets:
    ACT.checkPath(r['left'])
    ACT.checkPath(r['right'])

dump(ACT)

ACT.script(ACT.title, "BA3P - SNP calling and annotation", "BA3P")
ACT.begin(timestamp=False)

ACT.scene(1, "General configuration")
ACT.reportf("""Sample name: <b>{}</b><br>
Reference genome: <b>{}</b><br>
Bowtie2 index: <B>{}</b><br>
""".format(ACT.title, ACT.reference, ACT.btidx))
if ACT.mapfile != None:
    ACT.reportf("Chromosome mapping: <b>{}</b><br>".format(ACT.mapfile))
if ACT.snpeffdb != None:
    ACT.reportf("snpEff database: <b>{}</b><br>".format(ACT.snpeffdb))
ACT.reportf("Samples and input files:<br>")
ACT.table([ [r['name'], r['left'], r['right'] ] for r in SC.readsets],
          header=["Sample", "Left reads", "Right reads"],
          align="HLL")

#
# Before we start, check if reference file is indexed
#

fai = ACT.reference + ".fai" 
refdict = ACT.setFileExt(ACT.reference, ".dict")

if ACT.missingOrStale(fai, ACT.reference):
    print "Reference index {} does not exist or is out of date, creating it.".format(fai)
    ACT.shell("module load samtools; samtools faidx {}", ACT.reference)

if ACT.missingOrStale(refdict, ACT.reference):
    print "Reference dictionary {} does not exist, creating it.".format(refdict)
    ACT.submit("picard.qsub CreateSequenceDictionary R={} O={}".format(ACT.reference, refdict), done="dict.done")
    ACT.wait("dict.done")

# Ensure we don't have old .done files lying around
ACT.shell("rm -f *.done")

# And now run the pipeline....
sickleBowtie(ACT)
postBowtie(ACT)
ACT.doGATK()
ACT.postGATK()
ACT.countBAMs("bamcounts.csv", [(rs['name'], rs['bamfile']) for rs in ACT.sc.readsets], 'fixedReads')

ACT.scene(2, "Alignment")
ACT.reportf("""Reads were aligned to the reference genome using <b>bowtie2</b>, after trimming with <b>sickle</b>.
Aligned reads were processed with the <b>GATK</b> pipeline, to remove duplicates and fix mates.
The following table shows the original number of aligned reads and the number of reads after GATK processing for each sample.
Each sample name is linked to the corresponding BAM file.""")
ACT.table([ ["<a href='{}'>{}</a>".format(rs['bamfile'], rs['name']), rs['alignedReads'], rs['fixedReads'] ] for rs in ACT.sc.readsets], 
          header=['Sample', 'Aligned reads', 'Post-GATK reads'], align="HRR")

# SNP calling
if ACT.singleVCF:
    doVCFsingle(ACT)
else:
    doVCFmulti(ACT)

ACT.scene(3, "SNP calling")
if ACT.singleVCF:
    ACT.reportf("""SNPs were called using <b>freebayes</b> in single-VCF mode (all BAM files were combined
together for this step. A total of {} SNPs were identified. The VCF file can be downloaded using the link below.""", ACT.vcfsnps)
    ACT.file(ACT.singleVCF, description="VCF file containing all identified SNPs.")
else:
    ACT.reportf("""SNPs were called using <b>freebayes</b>. The following table shows the number of SNPs
identified in each sample. Each sample name is linked to the corresponding VCF file.""")
    ACT.table([ ["<a href='{}'>{}</a>".format(rs['vcf'], rs['name']), rs['vcfsnps'] ] for rs in ACT.sc.readsets],
              header=['Sample', 'SNPs'], align='HR')

lastscene = 4

# Run snpEff, if requested
if ACT.snpeffdb:
    lastscene = 5
    if ACT.singleVCF:
        snpeffSingle(ACT)
    else:
        snpeffMulti(ACT)

    ACT.scene(4, "SNP functional annotation")
    if ACT.singleVCF:
        ACT.reportf("""SNPs were annotated <b>snpEff</b> using the <b>{}</b> database. The following table provides links to the annotated VCF file
and to the snpEff reports (in CSV and HTML formats).""".format(ACT.snpeffdb))
        ACT.table([[linkify(ACT.snpeffVCF), linkify(ACT.csvreport), linkify(ACT.htmlreport)]],
                  header=['VCF', 'CSV report', 'HTML report'],
                  align="LLL")
    else:
        ACT.reportf("""SNPs were annotated through <b>snpEff</b> using the <b>{}</b> database. The following table provides links to the annotated VCF file
and to the snpEff reports (in CSV and HTML formats) for each sample.""".format(ACT.snpeffdb))
        ACT.table([ [rs['name'], linkify(rs['snpeff']), linkify(rs['csvreport']), linkify(rs['htmlreport']) ] for rs in ACT.sc.readsets],
                  header=['Sample', 'VCF', 'CSV report', 'HTML report'],
                  align="HLLL")

# Merge VCF files
if ACT.singleVCF:
    doVCFmergeSingle(ACT)
else:
    doVCFmergeMulti(ACT)

data = readReport(ACT)

ACT.scene(lastscene, "Combined SNP files")
ACT.reportf("""The VCF files for all samples were combined into a set of multi-fasta files (one for each chromosome).
The following table provides a link to the FASTA file for each chromosome, with the number of SNPs it contains.""")
ACT.table(data, header=['Chromosome', 'FASTA', 'SNPs'], align="HLR")

print "Script terminated, press Enter to delete unnecessary BAM files."
a=raw_input()

toDelete = ""
for rs in ACT.sc.readsets:
    toDelete += " " + rs['picard']
    toDelete += " " + rs['dupmark']
    toDelete += " " + rs['realigned']
ACT.execute("rm" + toDelete)


