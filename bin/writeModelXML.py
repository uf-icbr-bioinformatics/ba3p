#!/usr/bin/env python

# (c) 2015 A. Riva, DiBiG, UF ICBR Bioinformatics
# ariva@ufl.edu

import sys
import math
import random
import os.path

# Configuration object

class Configuration():
    locationsfile = False
    outfile = False
    inverseDescs = False
    # For the -g option
    nlocs = 0
    ndescs = 0

    def __init__(self, args):
        next = ""

        for a in args:
            if next == '-i':
                self.inverseDescs = True
                next = ""
            elif next == 'g0':
                self.locationsfile = a
                next = 'g1'
            elif next == 'g1':
                self.nlocs = int(a)
                next = 'g2'
            elif next == 'g2':
                self.ndescs = int(a)
                return
            elif a == '-g':
                next = 'g0'
            elif a == '-i':
                next = a
            elif self.locationsfile:
                self.outfile = a
            else:
                self.locationsfile = a

# A small class to represent descriptors

class Descriptor():
    name = ""
    description = ""
    values = []
    matrix = []

    def __init__(self, name, desc=""):
        self.name = name
        self.description = desc
        self.values = []
        self.matrix = []

    def createMatrix(self):
        nlocs = len(self.values)
        for i in range(nlocs):
            for j in range(nlocs):
                if i != j:
                    self.matrix.append(math.log(self.values[i] / self.values[j], 2))
        return self.matrix
    
    def inverseDescriptor(self):
        inv = Descriptor("{}-inv".format(self.name))
        inv.values = list(self.values)
        inv.matrix = [ -x for x in self.matrix ]
        return inv
    
def parseLocationsTable(filename):
    """Reads a table containing locations and descriptors from file `filename'. The file
should contain location names in the first column, and descriptors in all successive columns.
The function returns a tuple of two elements: locations and descriptors."""
    locations = []
    descriptors = []
    ndescriptors = 0
    hdr = True                  # Do we need to read the header?
    with open(filename, "r") as f:
        for line in f:
            parsed = line.strip("\r\n").split("\t")
            if hdr:
                descnames = parsed[1:]
                for d in descnames:
                    descriptors.append(Descriptor(d))
                    ndescriptors += 1
                hdr = False
            else:
                locations.append(parsed[0])
                for d, v in zip(descriptors, parsed[1:]):
                    d.values.append(float(v))
    return (locations, descriptors)

def loadLocations(filename, col=0):
    """Load a list of locations from file `filename'. The file
is assumed to be tab-delimited, and locations are read from column `col'."""
    locations = []
    with open(filename, "r") as f:
        for line in f:
            parsed = line.strip("\r\n").split("\t")
            locations.append(parsed[col])
    return locations

def writeGeneral(out, locations):
    sys.stderr.write("Writing General section\n")
    out.write('  <!-- name locations that are grouped -->\n')
    out.write('  <generalDataType id="geography">\n')
    for loc in locations:
        out.write('    <state code="{}"/>\n'.format(loc))
    out.write('  </generalDataType>\n')

def makeCoefficients(nd):
    """Returns a list containing `nd' values equal to 1/nd. Takes
care to ensure they sum to 1."""
    res = []
    x = 1.0/nd
    sum = 0
    for i in range(nd-1):
        res.append(str(x))
        sum += x
    res.append(str(1-sum))
    return res

def writeSubstitutionModel(out, locations, descriptors):
    ndescriptors = len(descriptors)
    coefficients = makeCoefficients(ndescriptors)
    indicators = [ str(1 - (i & 1)) for i in range(ndescriptors) ]

    sys.stderr.write("Writing Substitution Model\n")
    out.write("""  <glmSubstitutionModel id="originModel">
    <dataType idref="geography"/>
    <rootFrequencies>
      <frequencyModel>
        <dataType idref="geography"/>
          <frequencies>
            <parameter dimension="{}" value="{}" />
	  </frequencies>
      </frequencyModel>
    </rootFrequencies>
    <glmModel id="glmModel" family="logLinear"  checkIdentifiability="false">
      <independentVariables>
        <parameter id="glmTestCoefficients" value="{}" />
        <indicator>
          <parameter id="coefTestIndicator" value="{}" />
        </indicator>
        <designMatrix id="testDesignMatrix">
""".format(len(locations), "", " ".join(coefficients), " ".join(indicators)))

    for d in descriptors:
        out.write("""
          <!-- predictor: {} -->
          <parameter id="{}" value="{}" />
        """.format(d.description, d.name, " ".join([str(x) for x in d.matrix])))

    out.write("""
        </designMatrix>
      </independentVariables>
    </glmModel>
  </glmSubstitutionModel>
""")

def makeRewards(n, i):
    d = []
    for j in range(n):
        if j == i:
            d.append("1.0")
        else:
            d.append("0.0")
    return d

def makeFromLine(n, i):
    d = []
    for row in range(n):
        for col in range(n):
            if row == i:
                if col == i:
                    d.append("0")
                else:
                    d.append("1")
            else:
                d.append("0")
    return d

def makeToLine(n, i):
    d = []
    for row in range(n):
        for col in range(n):
            if col == i:
                if row == i:
                    d.append("0")
                else:
                    d.append("1")
            else:
                d.append("0")
    return d

def makeFromToLine(n, i, j):
    d = []
    for row in range(n):
        for col in range(n):
            if row == i and col == j:
                d.append("1")
            else:
                d.append("0")
    return d

def writeMarkov(out, locations):
    sys.stderr.write("Writing Markov section\n")
    n = len(locations)
    out.write("""  <markovJumpsTreeLikelihood id="geoTreeLikelihood" useAmbiguities="true" useUniformization="true" numberOfSimulants="1" saveCompleteHistory="true">
    <patterns idref="geoPatterns"/>
    <treeModel idref="treeModel"/>
    <siteModel idref="geoSiteModel"/>
    <complexSubstitutionModel idref="originModel"/>
      <rewards>
""")
    for i in range(n):
        data = makeRewards(n, i)
        out.write("""      <parameter id="{}_R" value="{}"/>\n""".format(locations[i], " ".join(data)))
    out.write("""      </rewards>\n""")

    # Now write ..._from lines
    for i in range(n):
        data = makeFromLine(n, i)
        out.write("""      <parameter id="{}_from" value="{}" />\n""".format(locations[i], " ".join(data)))

    # Now write ..._to lines
    for i in range(n):
        data = makeToLine(n, i)
        out.write("""      <parameter id="{}_to" value="{}" />\n""".format(locations[i], " ".join(data)))

    # Now write From_..._To_... lines
    for i in range(n):
        for j in range(i + 1, n):
            data = makeFromToLine(n, i, j)
            out.write("""      <parameter id="From_{}_To_{}" value="{}" />\n""".format(locations[i], locations[j], " ".join(data)))

    # Now write the 'reverse' From_..._To_... lines
    for i in range(n):
        for j in range(i + 1, n):
                data = makeFromToLine(n, j, i)
                out.write("""      <parameter id="From_{}_To_{}" value="{}" />\n""".format(locations[j], locations[i], " ".join(data)))

#def makeRandomDescriptors(ndescriptors, nlocations):
#    x = nlocations * (nlocations - 1)
#    return [ Descriptor("desc{}".format(i+1), "descriptor #{}".format(i+1), x) for i in range (ndescriptors) ]

def writeExampleFile(filename, nlocs, ndescs):
    with open(filename, "w") as out:
        out.write("Locations")
        for i in range(ndescs):
            out.write("\tdesc{}".format(i+1))
        out.write("\n")
        for l in range(nlocs):
            out.write("LOC{}".format(l+1))
            for i in range(ndescs):
                out.write("\t")
            out.write("\n")

def main(args):
    CONF = Configuration(args)
    out = sys.stdout

    if CONF.nlocs > 0:           # -g mode?
        writeExampleFile(CONF.locationsfile, CONF.nlocs, CONF.ndescs)
        return

    (locations, descriptors) = parseLocationsTable(CONF.locationsfile)
    for d in descriptors:
        d.createMatrix()
    if CONF.inverseDescs:
        descs2 = []
        for d in descriptors:
            descs2.append(d)
            d2 = d.inverseDescriptor()
            descs2.append(d2)
        descriptors = descs2

    nlocations = len(locations)
    ndescriptors = len(descriptors)
    sys.stderr.write("{} locations, {} descriptors read from `{}'.\n".format(nlocations, ndescriptors, CONF.locationsfile))
    if ndescriptors > 0:
        sys.stderr.write("Descriptors:\n")
        for d in descriptors:
            sys.stderr.write("  " + d.name + "\n")

    if CONF.outfile:
        out = open(CONF.outfile, "w")
    try:
        writeGeneral(out, locations)
        out.write("\n")
        writeSubstitutionModel(out, locations, descriptors)
        out.write("\n")
        writeMarkov(out, locations)
    finally:
        if CONF.outfile:
            sys.stderr.write("Output written to `{}'.\n".format(CONF.outfile))
            out.close()
    
if __name__ == "__main__":
    if len(sys.argv) >= 2:
        main(sys.argv[1:])
    else:
        progname = os.path.split(sys.argv[0])[1]
        sys.stderr.write("""{} - Write the linear model section of a BEAST XML file.

Usage: {} [-i] infile [outfile]
          [-g] outfile nlocs ndescs

The input file `infile' should contain location names, one per line. Additional columns
represent descriptors. The first line of the input file should contain descriptor names.

XML output will be written to `outfile' if specified, or to standard output.

The -i option causes an inverse descriptor to be added for each descriptors specified
in the input file.

If -g is specified, the program will create an empty locations file `outfile' for `nlocs'
locations and `ndescs' descriptors.\n""".format(progname, progname))
        sys.exit(1)

