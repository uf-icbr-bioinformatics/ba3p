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
                   'left': ACT.getConf("left", name),
                   'right': ACT.getConf("right", name)}

            if run['left'] == None:
                message("Configuration error: `left' undefined for sample `{}'".format(name))
                valid = False
            if run['right'] == None:
                message("Configuration error: `right' undefined for sample `{}'".format(name))
                valid = False

            self.runs.append(run)
            self.nruns += 1

def dumpmsc(msc):
    global message
    message("Title: {}", msc.title)
    message("{} samples:", msc.nruns)
    for run in msc.runs:
        message("  Name: {}", run['name'])
        message("  Left: {}", run['left'])
        message("  Right: {}", run['right'])
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
ACT.table([ [r['name'], r['left'], r['right'] ] for r in MSC.runs],
          header=["Sample", "Left reads", "Right reads"],
          align="HLL")

# Ensure we don't have old .done files lying around
ACT.shell("rm -f *.done")

#
# First of all run FastQC and sickle on fastq files
#

for run in MSC.runs:
    name = run['name']
    left = fixPath(run['left'])
    right = fixPath(run['right'])

    in1 = os.path.basename(left)
    in2 = os.path.basename(right)
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
    this = ACT.submit("sickle.qsub {} {} {} {}".format(left, right, run['sickle1'], run['sickle2']), done="sickle.done")
    fqc1 = ACT.submit("fastqc.qsub {} {}".format(left, run['fqcdir1']), done="fqc.done")
    fqc2 = ACT.submit("fastqc.qsub {} {}".format(right, run['fqcdir2']), done="fqc.done")
ACT.wait([('sickle.done', MSC.nruns), ('fqc.done', MSC.nruns * 2)])

for run in MSC.runs:
    fqc3 = ACT.submit("fastqc.qsub {} {}".format(run['sickle1'], run['fqcdir1a']), done="fqc2.done")
    fqc4 = ACT.submit("fastqc.qsub {} {}".format(run['sickle2'], run['fqcdir2a']), done="fqc2.done")
ACT.wait([('fqc2.done', MSC.nruns * 2)])

## Run spades on each sample

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
    infasta = smp['spadesdir'] + "contigs.fasta"
    outfasta = smp['name'] + ".spades.fasta"
    ACT.shell("cp {} Contigs/{}".format(infasta, outfasta))

## Run Mauve contig orderer

ACT.submit("mauve-contig-mover.qsub {}".format(ACT.reference))
