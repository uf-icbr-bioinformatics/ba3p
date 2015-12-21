#!/usr/bin/env python

# (c) 2015 A. Riva, DiBiG, UF ICBR Bioinformatics
# ariva@ufl.edu

import sys
import random
import os.path

# Configuration object

class Configuration():
    locationsfile = False
    outfile = False
    ndescriptors = 1

    def __init__(self, args):
        next = ""

        for a in args:
            if next == '-nd':
                self.ndescriptors = int(a)
                next = ""
            elif a == '-nd':
                next = a
            elif self.locationsfile:
                self.outfile = a
            else:
                self.locationsfile = a

# A small class to represent descriptors

class Descriptor():
    name = ""
    description = ""
    nvalues = 0
    values = []

    def __init__(self, name, desc="", nvalues=0):
        self.name = name
        self.description = desc
        self.nvalues = nvalues
        if nvalues > 0:
            self.values = [ str(random.uniform(-3.0, 3.0)) for i in range(nvalues) ]

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

def writeSubstitutionModel(out, locations, descriptors):
    ndescriptors = len(descriptors)
    coefficients = ["0.1" for i in range(ndescriptors) ]
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
""".format(d.description, d.name, " ".join(d.values)))

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
    

def makeRandomDescriptors(ndescriptors, nlocations):
    x = nlocations * (nlocations - 1)
    return [ Descriptor("desc{}".format(i+1), "descriptor #{}".format(i+1), x) for i in range (ndescriptors) ]
        
def main(args):
    CONF = Configuration(args)
    out = sys.stdout

    locations = loadLocations(CONF.locationsfile)
    nlocations = len(locations)
    sys.stderr.write("{} locations read from `{}'.\n".format(nlocations, CONF.locationsfile))

    if CONF.outfile:
        out = open(CONF.outfile, "w")
    try:
        writeGeneral(out, locations)
        out.write("\n")
        writeSubstitutionModel(out, locations, makeRandomDescriptors(CONF.ndescriptors, nlocations))
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

Usage: {} [-nd N] infile [outfile]

The input file `infile' should contain location names, one per line. XML output
will be written to `outfile' if specified, or to standard output.

Use the -nd argument to specify the number of descriptors (1 by default).\n""".format(progname, progname))
        sys.exit(1)

