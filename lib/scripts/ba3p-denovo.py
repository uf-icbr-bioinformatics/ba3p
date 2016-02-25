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

class ba3pVOPT():
    valid = True                # Set to False in case of configuration errors
    title = "denovo"            # Title of run
    runs = []                   # List of runs
    nruns = 0                   # Number of runs

    def __init__(self, ACT):
        conffile = ACT.Arguments[0]
        conf = ACT.loadConfiguration(conffile)
        self.title = ACT.getConf("title")

        runnames = ACT.getConf("samples").split(",")
        runnames = [ s.strip() for s in runnames ]

        for name in runnames:
            run = {'name': name,
                   'fastq1': ACT.getConf("fastq1", name),
                   'fastq2': ACT.getConf("fastq2", name)}

            if run['fastq1'] == None:
                myutils.message("Configuration error: `fastq1' undefined for sample `{}'".format(name))
                valid = False
            if run['fastq2'] == None:
                myutils.message("Configuration error: `fastq2' undefined for sample `{}'".format(name))
                valid = False

            self.runs.append(run)
            self.nruns += 1

def dumpmsc(msc):
    global message
    message("Title: {}", msc.title)
    message("{} samples:", msc.nruns)
    for run in msc.runs:
        message("  Name: {}", run['name'])
        message("  Fastq1: {}", run['fastq1'])
        message("  Fastq2: {}", run['fastq2'])
        message("")

MSC = ba3pVOPT(ACT)

if not MSC.valid:
    message("Script cannot be executed due to configuration errors.")
    exit
dumpmsc(MSC)

# Script definition 

ACT.script(MSC.title, "BA3P - De-novo Assembly", "BA3P")
ACT.begin(timestamp=False)

ACT.scene(1, "General configuration")
ACT.reportf("""Experiment name: <b>{}</b><br>
""".format(MSC.title))
ACT.reportf("Samples and input files:<br>")
ACT.table([ [r['name'], r['fastq1'], r['fastq2'] ] for r in MSC.runs],
          header=["Sample", "Left reads", "Right reads"],
          align="HLL")

# Ensure we don't have old .done files lying around
ACT.shell("rm -f *.done")

#
# First of all run FastQC and sickle on fastq files
#

for run in MSC.runs:
    name = run['name']
    fastq1 = fixPath(run['fastq1'])
    fastq2 = fixPath(run['fastq2'])

    in1 = os.path.basename(fastq1)
    in2 = os.path.basename(fastq2)
    run['fqcdir1'] = ACT.setFileExt(in1, ".before.fqc", remove=[".fastq", ".gz"])
    run['fqcdir2'] = ACT.setFileExt(in2, ".before.fqc", remove=[".fastq", ".gz"])
    run['fqcdir1a'] = ACT.setFileExt(in1, ".after.fqc", remove=[".fastq", ".gz"])
    run['fqcdir2a'] = ACT.setFileExt(in2, ".after.fqc", remove=[".fastq", ".gz"])
    ACT.mkdir(run['fqcdir1'])
    ACT.mkdir(run['fqcdir2'])
    ACT.mkdir(run['fqcdir1a'])
    ACT.mkdir(run['fqcdir2a'])
    run['sickle1'] = ACT.setFileExt(in1, ".sickle.fastq", remove=[".fastq", ".gz"])
    run['sickle2'] = ACT.setFileExt(in2, ".sickle.fastq", remove=[".fastq", ".gz"])
    this = ACT.submit("sickle.qsub {} {} {} {}".format(fastq1, fastq2, run['sickle1'], run['sickle2']), done="sickle.done")
    fqc1 = ACT.submit("fastqc.qsub {} {}".format(fastq1, run['fqcdir1']), done="fqc.done")
    fqc2 = ACT.submit("fastqc.qsub {} {}".format(fastq2, run['fqcdir2']), done="fqc.done")
ACT.wait([('sickle.done', MSC.nruns), ('fqc.done', MSC.nruns * 2)])

for run in MSC.runs:
    fqc3 = ACT.submit("fastqc.qsub {} {}".format(run['sickle1'], run['fqcdir1a']), done="fqc2.done")
    fqc4 = ACT.submit("fastqc.qsub {} {}".format(run['sickle2'], run['fqcdir2a']), done="fqc2.done")
ACT.wait([('fqc2.done', MSC.nruns * 2)])

## Run spads on each sample

nspades = 0
for smp in MSC.samples:
    outdir = smp['name'] + ".spades/"
    smp['spadesdir'] = outdir
    ACT.mkdir(outdir)
    ACT.submit("spades.qsub {} {} {}".format(outdir, smp['sickle1'], smp['sickle2']), done="spades.@.done")
    nspades = nspades + 1
ACT.wait(("spades.@.done", nspades))

## Create contigs directory and copy fasta files into it

ACT.mkdir("Contigs")
for smp in MSC.samples:
    fasta = smp['spadesdir'] + "..."
    ACT.shell("cp {} Contigs/".format(fasta))

## Run Mauve contig orderer

ACT.submit("mauve-contig-mover.qsub {}".format(ACT.reference))
