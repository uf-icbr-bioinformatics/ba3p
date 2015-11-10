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
    fqcdir1 = in1 + ".before.fqc/"
    fqcdir2 = in2 + ".before.fqc/"
    ACT.mkdir(fqcdir1)
    ACT.mkdir(fqcdir2)
    run['sickle1'] = ACT.setFileExt(in1, ".sickle.fastq", remove=[".fastq", ".gz"])
    run['sickle2'] = ACT.setFileExt(in2, ".sickle.fastq", remove=[".fastq", ".gz"])
    this = ACT.submit("sickle.qsub {} {} {} {}".format(fastq1, fastq2, run['sickle1'], run['sickle2']), done="sickle.done")
    fqc1 = ACT.submit("fastqc.qsub {} {}".format(fastq1, fqcdir1), done="fqc1.done")
    fqc2 = ACT.submit("fastqc.qsub {} {}".format(fastq2, fqcdir2), done="fqc2.done")
     # print "Submitted: " + this
ACT.wait(('sickle.done', MSC.nruns), ('fqc1.done', MSC.nruns), ('fcq2.done', MSC.runs))

#
# Once all the bowties are done, we can collect stats (number of aligned reads),
# and then proceed with the GATK pipeline, followed by SNP calling.
#

#for run in MSC.runs:
#    ACT.submit("gzip.qsub {}".format(run['sickle1']))
#    ACT.submit("gzip.qsub {}".format(run['sickle2']))

