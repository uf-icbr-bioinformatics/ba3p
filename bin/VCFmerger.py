#!/usr/bin/env python

#################################
#
# (c) 2015, A. Riva
# DiBiG - ICBR Bioinformatics
#
#################################

import sys
import os.path

# Globals

QUALFIELD = "standard_min_confidence_threshold_for_calling="

IMPACTS = {"HIGH": 4, "MODERATE": 3, "LOW": 2, "MODIFIER": 1}

ANNOTATIONS = [
    ("HIGH", "chromosome_number_variation", ""),
    ("HIGH", "chromosome", "A large parte (over 1%) of the chromosome was deleted."),
    ("HIGH", "exon_loss_variant", "A deletion removes the whole exon."),
    ("HIGH", "frameshift_variant", "Insertion or deletion causes a frame shift e.g.: An indel size is not multple of 3"),
    ("HIGH", "rare_amino_acid_variant", "The variant hits a rare amino acid thus is likely to produce protein loss of function"),
    ("HIGH", "splice_acceptor_variant", "The variant hits a splice acceptor site (defined as two bases before exon start, except for the first exon)."),
    ("HIGH", "splice_donor_variant", "The variant hits a Splice donor site (defined as two bases after coding exon end, except for the last exon)."),
    ("HIGH", "start_lost", "Variant causes start codon to be mutated into a non-start codon. e.g.: aTg/aGg, M/R"),
    ("HIGH", "stop_gained", "Variant causes a STOP codon e.g.: Cag/Tag, Q/*"),
    ("HIGH", "stop_lost", "Variant causes stop codon to be mutated into a non-stop codon e.g.: Tga/Cga, */R"),
    ("HIGH", "transcript_ablation", ""),

    ("MODERATE", "3_prime_UTR_truncation+exon_loss", "The variant deletes an exon which is in the 3'UTR of the transcript"),
    ("MODERATE", "5_prime_UTR_truncation+exon_loss_variant", "The variant deletes an exon which is in the 5'UTR of the transcript"),
    ("MODERATE", "coding_sequence_variant", "One or many codons are changed e.g.: An MNP of size multiple of 3"),
    ("MODERATE", "disruptive_inframe_deletion", "One codon is changed and one or more codons are deleted e.g.: A deletion of size multiple of three, not at codon boundary"),
    ("MODERATE", "disruptive_inframe_insertion", "One codon is changed and one or many codons are inserted e.g.: An insert of size multiple of three, not at codon boundary"),
    ("MODERATE", "inframe_deletion", "One or many codons are deleted e.g.: A deletion multiple of three at codon boundary"),
    ("MODERATE", "inframe_insertion", "One or many codons are inserted e.g.: An insert multiple of three in a codon boundary"),
    ("MODERATE", "missense_variant", "Variant causes a codon that produces a different amino acid e.g.: Tgg/Cgg, W/R"),
    ("MODERATE", "sequence_feature+exon_loss_variant", "A 'NextProt' based annotation. Details are provided in the 'feature type' sub-field (ANN), or in the effect details (EFF)."),
    ("MODERATE", "splice_region_variant", "A variant affecting putative (Lariat) branch point from U12 splicing machinery, located in the intron."),
    ("MODERATE", "regulatory_region_ablation", ""),
    ("MODERATE", "TFBS_ablation", ""),

    ("LOW", "5_prime_UTR_premature_start_codon_gain_variant", "A variant in 5'UTR region produces a three base sequence that can be a START codon."),
    ("LOW", "initiator_codon_variant", "Variant causes start codon to be mutated into another start codon (the new codon produces a different AA). e.g.: Atg/Ctg, M/L (ATG and CTG can be START codons)"),
    ("LOW", "splice_region_variant", "A sequence variant in which a change has occurred within the region of the splice site, either within 1-3 bases of the exon or 3-8 bases of the intron."),
    ("LOW", "start_retained", "Variant causes start codon to be mutated into another start codon. e.g.: Ttg/Ctg, L/L (TTG and CTG can be START codons)"),
    ("LOW", "stop_retained_variant", "Variant causes stop codon to be mutated into another stop codon. e.g.: taA/taG, */*"),
    ("LOW", "synonymous_variant", "Variant causes a codon that produces the same amino acid. e.g.: Ttg/Ctg, L/L"),

    ("MODIFIER", "3_prime_UTR_variant", "Variant hits 3'UTR region"),
    ("MODIFIER", "5_prime_UTR_variant", "Variant hits 5'UTR region"),
    ("MODIFIER", "coding_sequence_variant", "The variant hits a CDS."),
    ("MODIFIER", "conserved_intergenic_variant", "The variant is in a highly conserved intergenic region"),
    ("MODIFIER", "conserved_intron_variant", "The variant is in a highly conserved intronic region"),
    ("MODIFIER", "downstream_gene_variant", "Downstream of a gene (default length: 5K bases)"),
    ("MODIFIER", "exon_variant", "The variant hits an exon (from a non-coding transcript) or a retained intron."),
    ("MODIFIER", "gene_variant", "The variant hits a gene."),
    ("MODIFIER", "feature_elongation", ""),
    ("MODIFIER", "feature_truncation", ""),
    ("MODIFIER", "intergenic_region", "The variant is in an intergenic region"),
    ("MODIFIER", "intragenic_variant", "The variant hits a gene, but no transcripts within the gene"),
    ("MODIFIER", "intron_variant", "Variant hits and intron. Technically, hits no exon in the transcript."),
    ("MODIFIER", "mature_miRNA_variant", ""),
    ("MODIFIER", "miRNA", "Variant affects a miRNA"),
    ("MODIFIER", "NMD_transcript_variant", ""),
    ("MODIFIER", "non_coding_transcript_exon_variant", ""),
    ("MODIFIER", "non_coding_transcript_variant", ""),
    ("MODIFIER", "regulatory_region_amplification", ""),
    ("MODIFIER", "regulatory_region_variant", "The variant hits a known regulatory feature (non-coding)."),
    ("MODIFIER", "TF_binding_site_variant", ""),
    ("MODIFIER", "TFBS_amplification", ""),
    ("MODIFIER", "transcript_amplification", ""),
    ("MODIFIER", "transcript_variant", "The variant hits a transcript."),
    ("MODIFIER", "upstream_gene_variant", "Upstream of a gene (default length: 5K bases)")
]

# Utils

def safeReadfloat(s):
    try:
        return float(s)
    except ValueError:
        return None

def safeReadInt(s):
    try:
        return int(s)
    except ValueError:
        return None

def parsePair(s):
    p = s.find("/")
    if p == -1:
        return ()
    else:
        a = safeReadInt(s[0:p])
        b = safeReadInt(s[p+1:])
        return (a, b)

# Main classes

class VCFSNP:
    position = 0
    nseen = 0
    reference = ""
    alternate = ""
    genotypes = []
    snpeff = None

    def __init__(self, position, nfiles):
        self.position = position
        self.alleles   = [ "-" for i in range(nfiles) ]
        self.genotypes = [ 0 for i in range(nfiles) ]

    def findHighestImpact(self):
        nimp = 0

        if self.snpeff != None:
            for ann in self.snpeff.annots:
                # print ann.impact
                thisn = IMPACTS[ann.impact]
                if thisn > nimp:
                    nimp = thisn
        return nimp

    def findHighestImpactString(self):
        imp = ""
        nimp = 0

        if self.snpeff != None:
            for ann in self.snpeff.annots:
                # print ann.impact
                thisn = IMPACTS[ann.impact]
                if thisn > nimp:
                    imp = ann.impact
                    nimp = thisn
        return imp

    def impactLessThan(self, wanted):
        """Returns True if the highest impact of this SNP is less than `wanted'."""
        nimp = self.findHighestImpact()
        # print "Comparing {} with {}.".format(nimp, wanted)
        return (nimp < wanted)

class VCFfile:
    filename = None
    samples = []
    nsamples = 0
    qual = None

    def __init__(self, filename):
        self.filename = filename
        if not os.path.isfile(filename):
            sys.stderr.write("File {} does not exist or is not readable.\n".format(filename))
            return None
        with open(filename, "r") as f:
            for line in f:
                if len(line) > 0:                                         # skip empty lines (there should not be any)
                    if line[0:2] == "##":                                 # still in header?
                        self.parseQualFromVCF(line)
                    elif line[0] == "#":
                        line = line.rstrip("\n").split("\t")
                        self.samples = line[9:]
                        self.nsamples = len(self.samples)
                    else:
                        return
        
    def parseQualFromVCF(self, line):
        """If QUALFIELD is present in `line', extract the quality value that
follows it and set self.qual to it, signaling that quality was retrieved from
the VCF file."""
        p1 = line.find(QUALFIELD)
        if p1 >= 0:
            p2 = line.find(" ", p1)
            if p2 >= 0:
                self.qual = safeReadfloat(line[p1+len(QUALFIELD):p2]) # this will return None in case of parsing errors

class VCFparser:
    chromosomes = []            # list of chromosomes in order seen
    snps = {}                   # dictionary mapping chromosomes to lists of VCFSNPs
    filenames = []
    files = []
    nfiles = 0
    nsamples = 0                # Total number of samples
    qual = 30
    qualHow = "d"               # quality set by default

    def __init__(self, filenames, qual=None):
        self.filenames = filenames
        self.nfiles = len(filenames)
        self.referenceBases = {}
        self.snpSets = []
        if qual != None:
            self.qual = qual
            self.qualHow = "a"  # quality was set through command-line argument
        self.files = [ VCFfile(f) for f in filenames ]
        for vf in self.files:
            self.nsamples += vf.nsamples

    def parseQualFromVCF(self, line):
        """If QUALFIELD is present in `line', extract the quality value that
follows it and set self.qual to it, signaling that quality was retrieved from
the VCF file."""
        p1 = line.find(QUALFIELD)
        if p1 >= 0:
            p2 = line.find(" ", p1)
            if p2 >= 0:
                q = safeReadfloat(line[p1+len(QUALFIELD):p2])
                if q != None:
                    self.qual = q
                    self.qualHow = "v"

    def extractSnpEff(self, line):
        parsed = line.split(";")
        for p in parsed:
            if p[0:4] == "ANN=":
                return snpEffField(p[4:])
        return None

    def parseAllVCF(self):
        "Parse all VCF files."
        idx = 0
        for f in self.files:
            print "Parsing {}...".format(f.filename)
            self.parseOneVCF(f, idx)
            idx += f.nsamples

    def parseOneVCF(self, vf, idx):
        """Parse SNP information from VCFfile `vf', storing its genotypes starting at position `idx' in the output vector."""
        inHeader = True
        thisChrom = ""          # Current chromosome
        thisSNPs = None         # List of SNPs for current chromosome
        with open(vf.filename, "r") as f:
            for line in f:
                if len(line) > 0:                                             # skip empty lines (there should not be any)
                    if inHeader:
                        if line[0:2] == "##":                                 # still in header?
                            if (self.qualHow == "d"):
                                self.parseQualFromVCF(line)
                        else:
                            inHeader = False                                  # header done (this also skips the column headers)
                            # print "End of header."
                            how = False
                            if self.qualHow == "d":
                                how = "(default)"
                                self.qualHow = "D"
                            elif self.qualHow == "a":
                                how = "(command-line)"
                                self.qualHow = "A"
                            elif self.qualHow == "v":
                                how = "(VCF header)"
                                self.qualHow = "V"
                            if how:
                                print "Quality score threshold: {} {}".format(self.qual, how)
                            
                    else:       # Data line
                        line = line.split("\t")
                        if len(line) > 6: # otherwise something's really wrong with this VCF
                            pfield = line[6]
                            qfield = float(line[5])
                            if ((pfield == "PASS") or (pfield == ".")) and (qfield >= self.qual):
                                chrom = line[0]
                                pos = int(line[1])
                                ref = line[3].strip()
                                snp = line[4].strip()
                                alldata  = line[9:] # Columns containing genotypes

                                if len(ref) == 1 and len(snp) == 1:              # only consider SNPs
                                    if not chrom == thisChrom:                   # new chromosome?
                                        if not chrom in self.chromosomes:
                                            self.chromosomes.append(chrom)       # add this chromosome to the list, if not already there
                                        if chrom in self.snps:                   # have we seen this chrom before (in a previous file)?
                                            thisSNPs = self.snps[chrom]
                                        else:
                                            thisSNPs = {}                        # create new list of SNPs
                                            self.snps[chrom] = thisSNPs          # and add it to the dictionary
                                        thisChrom = chrom
                                        # print "new chromosome {} found.".format(chrom)
                                    if pos in thisSNPs:                          # do we already have a SNP in this position?
                                        s = thisSNPs[pos]
                                    else:
                                        s = VCFSNP(pos, self.nsamples)
                                        # print "Adding {} to {}.".format(s, self.snps)
                                        s.reference = ref
                                        s.alternate = snp
                                        s.snpeff = self.extractSnpEff(line[7])
                                        self.snps[thisChrom][pos] = s


                                    for i in range(vf.nsamples):
                                        dat = alldata[i]
                                        if dat.find(":") > 0:
                                            s.nseen += 1
                                            gt = dat.split(":")[0]
                                            gt1 = int(gt[0])
                                            gt2 = int(gt[2])
                                            # print "{}: {} + {} = {}".format(pos, gt1, gt2, gt1+gt2)
                                            s.genotypes[idx + i] = gt1 + gt2 + 1    # 1=ref hom, 2=het, 3=alt hom

    def summarize(self):
        print "Chromosomes: {}".format(", ".join(self.chromosomes))
        print "SNPs:"
        for chrom in self.chromosomes:
            print "  {}  {}".format(chrom, len(self.snps[chrom]))

    def filterSNPs(self, missing=False, impact=False):
        """Remove SNPs that were not seen in all files (if `missing' is True) or that have a highest impact less than `impact'."""
        removedMissing = 0
        removedImpact = 0
        wanted = self.nsamples
        for chrom in self.chromosomes:
            new = {}
            for (pos, v) in self.snps[chrom].iteritems():
                keep = True
                if missing and v.nseen < wanted:
                    removedMissing += 1
                    keep = False
                if impact and v.impactLessThan(impact):
                    removedImpact += 1
                    keep = False
                if keep:
                    new[pos] = v
            self.snps[chrom] = new
        if removedMissing > 0:
            print "{} SNPs removed because of missing genotypes.".format(removedMissing)
        if removedImpact > 0:
            print "{} SNPs removed because of low impact.".format(removedImpact)

#
# Parsing of snpEff annotations
#

class snpEffAnnot():
    allele = ""
    effect = ""
    impact = ""
    geneName = ""
    geneID = ""
    featureType = ""
    featureID = ""
    biotype = ""
    rank = ()
    HGVS_c = ""
    HGVS_p = ""
    cDNA = ()
    CDS = ()
    protein = ""
    distances = ""
    info = ""

    def __init__(self, fields):
        self.allele = fields[0]
        self.effect = fields[1]
        self.impact = fields[2]
        self.geneName = fields[3]
        self.geneID = fields[4]
        self.featureType = fields[5]
        self.featureID = fields[6]
        self.biotype = fields[7]
        self.rank = parsePair(fields[8])
        self.HGVS_c = fields[9]
        self.HGVS_p = fields[10]
        self.cDNA = parsePair(fields[11])
        self.CDS = parsePair(fields[12])
        self.protein = fields[13]
        self.distances = fields[14]
        self.info = fields[15]

class snpEffField():
    annots = []

    def __init__(self, field):
        self.annots = []
        entries = field.split(",")
        for e in entries:
            pe = e.split("|")
            ann = snpEffAnnot(pe)
            self.annots.append(ann)

class VCFmerger():
    outpattern = "out-{}.fa"
    bothAlleles = False
    missingChar = False         # If set to a character, missing SNPs will be indicated by this character instead of reference
    removeMissing = False       # If true, SNPs will be printed only if they have alleles in all samples
    quality = None
    impact = None               # Desired minimum impact
    nimpact = False             # Numeric value of desired minimum impact
    annotations = []
    snpsfile = None
    reportfile = None

    def showConf(self, nfiles):
        print """
#############################
# Input VCF files: {}
# outPattern: {}
# missingChar: {}
# bothAlleles: {}
# removeMissing: {}
# quality: {}
# impact: {}
# annotations: {}
#############################
""".format(nfiles, self.outpattern, self.missingChar, self.bothAlleles, self.removeMissing, self.quality, self.impact, ",".join(self.annotations))

    def parseOptions(self, args):
        next = ""
        filenames = []

        for a in args:
            if a in ["-m", "--missingChar"]:
                next = "-m"
            elif a in ["-o", "--outPattern"]:
                next = "-o"
            elif next == "-m":
                self.missingChar = a
                next = ""
            elif next == "-o":
                self.outpattern = a
                next = ""
            elif next == "-q":
                self.quality = safeReadfloat(a)
                next = ""
            elif next == "-i":
                if a in IMPACTS:
                    self.impact = a
                    self.nimpact = IMPACTS[a]
                next = ""
            elif next == "-a":
                self.annotations = a.split(",")
                next = ""
            elif next == "-s":
                self.snpsfile = a
                next = ""
            elif next == "-R":
                self.reportfile = a
                next = ""
            elif a in ["-b", "--bothAlleles"]:
                self.bothAlleles = True
                next = ""
            elif a in ["-r", "--removeMissing"]:
                self.removeMissing = True
                next = ""
            elif a in ["-q", "--quality"]:
                next = "-q"
            elif a in ["-i", "--impact"]:
                next = "-i"
            elif a in ["-a", "--annotation"]:
                next = "-a"
            elif a in ["-s", "--saveSNPs"]:
                next = "-s"
            elif a in ["-R", "--report"]:
                next = "-R"
            else:
                filenames.append(a)

        self.showConf(len(filenames))
        return filenames

    def writeFastaFiles(self, parser): # parser is a VCFparser object
        op = self.outpattern
        if op.find("{") >= 0:
            for chrom in parser.chromosomes:
                outfile = op.format(chrom)
                print "Writing output file {}...".format(outfile)
                with open(outfile, "w") as f:
                    self.writeOneChrom(parser, chrom, f)
        else:
            with open(op, "w") as f:
                for chrom in parser.chromosomes:
                    self.writeOneChrom(parser, chrom, f)

    def writeOneChrom(self, parser, chrom, f):
        idx = 0
        snplist = parser.snps[chrom]
        positions = sorted(snplist)
        # print positions
        for vf in parser.files:                             # Go through all files
            for i in range(vf.nsamples):                    # And all samples in each file. i is index of sample in file.
                f.write(">{}\n".format(vf.samples[idx + i]))      # Write sample name
                for p in positions:                         # Loop through all positions
                    v = snplist[p]
                    ref = v.reference
                    alt = v.alternate
                    gt = v.genotypes[idx + i]
                    # print v.genotypes
                    # print (ref, alt, gt)

                    if gt == 0:
                        if self.missingChar:
                            a = self.missingChar
                        else:
                            a = ref
                        if self.bothAlleles:
                            c = a + a
                        else:
                            c = a

                    elif gt == 1:
                        if self.bothAlleles:
                            c = ref + ref
                        else:
                            c = ref

                    elif gt == 2:
                        if self.bothAlleles:
                            c = ref + al
                        else:
                            c = alt

                    elif gt == 3:
                        if self.bothAlleles:
                            c = alt + alt
                        else:
                            c = alt

                    f.write(c)
                f.write("\n")
            idx += vf.nsamples

    def dumpSNPs(self, parser):
        filename = self.snpsfile
        with open(filename, "w") as out:
            out.write("Chrom\tPos\tRef\tAlt\tSeen\tImpact\n")
            for chrom in parser.chromosomes:
                for idx in range(parser.nfiles):
                    snplist = parser.snps[chrom]
                    positions = sorted(snplist)
                    for p in positions:
                        v = snplist[p]
                        ref = v.reference
                        alt = v.alternate
                        out.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom, p, ref, alt, v.nseen, v.findHighestImpactString()))

    def writeReport(self, parser):
        filename = self.reportfile
        with open(filename, "w") as out:
            out.write("Chrom\tFilename\tSNPs\n")
            for chrom in parser.chromosomes:
                outfile = self.outpattern.format(chrom)
                nsnps = len(parser.snps[chrom])
                out.write("{}\t{}\t{}\n".format(chrom, outfile, nsnps))
            
def main(args):
    m = VCFmerger()
    files = m.parseOptions(sys.argv[1:])
    p = VCFparser(files, qual=m.quality)
    p.parseAllVCF()

    #snps = p.snps['CFSAN001992']
    #print snps[141].genotypes
    #print snps[197].genotypes

    p.filterSNPs(missing=m.removeMissing, impact=m.nimpact)
    p.summarize()
    m.writeFastaFiles(p)
    if m.snpsfile != None:
        m.dumpSNPs(p)
    if m.reportfile != None:
        m.writeReport(p)

def usage():
    print """Usage: VCFmerger.py [options] vcffiles ...

Options:
  -m|--missingChar C         Use character C to represent unknown alleles (otherwise
                             reference allele will be used).
  -o|--outPattern P          Use string P as the name of the output file(s). If P contains
                             the sequence {}, it will be replaced by each chromosome
                             name, otherwise output will be written to a single file.
  -b|--bothAlleles           If specified, the program will print two bases for each
                             SNP, according to the genotype (AA, AB, or BB).
  -r|--removeMissing         If specified, the program will only print SNPs that 
                             appear in all the VCF files.
  -q|--quality Q             Use Q as the quality threshold. If not specified, the quality
                             threshold will be read from the VCF file, if present, otherwise
                             it will be set to the default value of 30.
  -R|--report R              Write the list of output files produced to a report file R in 
                             tab-delimited format. Each line contains three columns: 
                             chromosome, filename, number of SNPs.
  -i|--impact I              If snpEff annotations are present, only retain SNPs with a
                             putative impact I or higher. I should be one of HIGH, MODERATE,
                             LOW, or MODIFIER.
  -a|--annotation A[,A...]   If snpEff annotations are present, only retain SNPs having an
                             annotation matching one or more of the supplied annotations A.
                             See http://snpeff.sourceforge.net/VCFannotationformat_v1.0.pdf
                             for a list of all valid annotations, or use the -l option.
  -s|--saveSNPs S            Save SNPs to tab-delimited file S.
  -l|--listEff               List all valid values for snpEff impact and annotation.
                             When the long form is used, also prints a short description
                             of each annotation (if available).

"""

def listAnnotations(full=False):
    sys.stdout.write("""Impact   Annotation
-------- ---------
""")
    for ann in ANNOTATIONS:
        sys.stdout.write("{:8} {}\n".format(*ann))
        if full and ann[2] != "":
            sys.stdout.write("         {}\n".format(ann[2]))
    
if __name__ == "__main__":
    args = sys.argv[1:]
    if len(args) == 0 or args[0] in ["-h", "--help"]:
        usage()
    elif args[0] in ["-l", "--listEff"]:
        listAnnotations(full=(args[0]=="--listEff"))
    else:
        main(args)

