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
    e.write("""sam-remap.py - converts chromosome names in a SAM file.

Usage: sam-remap.py mapfile

Details:
  mapfile should be a tab-delimited file with two columns. The first column
contains chromosome names as they currently appear in the SAM file, and the
second column contains the new chromosome names.

This command is meant to work as a filter: it reads a SAM file from standard
input and writes a new SAM file to standard output. For example:

  > cat existing.sam | sam-remap.py map.txt > converted.sam

To convert chromosome names in a BAM file use something like the following:

  > samtools view -h existing.bam | sam-remap.py map.txt > temp.sam
  > samtools view -b -o converted.bam temp.sam; rm temp.sam

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
        parsed = line.split("\t")
        if line[0] == '@':
            # Header line
            if parsed[0] == "@SQ":
                orig = parsed[1][3:]
                if orig in trmap:
                    parsed[1] = "SN:" + trmap[orig]
        else:
            orig = parsed[2]
            if orig in trmap:
                parsed[2] = trmap[orig]
            orig = parsed[6]
            if orig in trmap:
                parsed[6] = trmap[orig]
            
        sys.stdout.write("\t".join(parsed))



    
