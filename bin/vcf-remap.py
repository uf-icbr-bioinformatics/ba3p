#!/usr/bin/env python

import sys

def parseLine(line):
    return line.rstrip("\n\r").split("\t")

def readMap(filename):
    trmap = {}
    with open(filename, "r") as f:
        for line in f:
            parsed = parseLine(line)
            trmap[parsed[0]] = parsed[1]
    return trmap

def usage():
    e = sys.stderr
    e.write("""vcf-remap.py - converts chromosome names in a VCF file.

Usage: vcf-remap.py mapfile

Details:
  mapfile should be a tab-delimited file with two columns. The first column
contains chromosome names as they currently appear in the VCF file, and the
second column contains the new chromosome names.

This command is meant to work as a filter: it reads a VCF file from standard
input and writes a new VCF file to standard output. For example:

  > cat existing.vcf | vcf-remap.py map.txt > converted.vcf

""")

if __name__=="__main__":

    if len(sys.argv) != 2:
        usage()
        exit(1)

    mapfile = sys.argv[1]
    trmap = readMap(mapfile)
    # print "map: {}".format(trmap)

    while True:
        line = sys.stdin.readline()
        if not line:
            exit(0)
        if line[0] == '#':
            sys.stdout.write(line)
        else:
            parsed = line.split("\t")
            orig = parsed[0]
            if orig in trmap:
                parsed[0] = trmap[orig]
            sys.stdout.write("\t".join(parsed))



    
