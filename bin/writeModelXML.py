#!/usr/bin/env python

# (c) 2015 A. Riva, DiBiG, UF ICBR Bioinformatics
# ariva@ufl.edu

import sys
import os.path

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
    

def main(args):
    outfile = False
    out = sys.stdout
    locationsfile = args[0]

    if len(args) > 1:
        outfile = args[1]

    locations = loadLocations(locationsfile)
    sys.stderr.write("{} locations read from `{}'.\n".format(len(locations), locationsfile))

    if outfile:
        out = open(outfile, "w")
    try:
        writeGeneral(out, locations)
        out.write("\n")
        writeMarkov(out, locations)
    finally:
        if outfile:
            sys.stderr.write("Output written to `{}'.\n".format(outfile))
            out.close()
    
if __name__ == "__main__":
    if len(sys.argv) >= 2:
        main(sys.argv[1:])
    else:
        sys.stderr.write("""{} - Write the linear model section of a BEAST XML file.

Usage: {} infile [outfile]

The input file `infile' should contain location names, one per line. XML output
will be written to `outfile' if specified, or to standard output.
""".format(sys.argv[0], sys.argv[0]))
        sys.exit(1)

