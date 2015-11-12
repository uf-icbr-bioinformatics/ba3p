#!/usr/bin/env python

# SWIPT - Sliding WIndow Phi Test
# (c) 2015 A. Riva, DiBiG, ICBR Bioinformatics, University of Florida
#          and M. Salemi, University of Florida
#

import os
import sys
import tempfile
import subprocess
from Bio import SeqIO

# Globals

class Globals():
    PHI = "/apps/dibig_ba3p/1.0/bin/Phi"   # -p
    PHIPERM = False                        # -phip
    PHIEXTRA = False                       # -phio
    WINSIZE = 200                          # -w
    WINSTEP = 40                           # -s
    PVALTHR = False                        # -t
    INFILE = False
    OUTFILE = False                        # -o

# Utils

def safeReadfloat(s):
    try:
        return float(s)
    except ValueError:
        return None

def tmpfile():
    (fd, path) = tempfile.mkstemp(dir=".")
    os.close(fd)
    return path

# Application functions

def invokePhi(cmdline):
    try:
        return subprocess.check_output(cmdline)
    except subprocess.CalledProcessError as e:
        return e.output
    except:
        return False

def findFloatTarget(lines, target):
    p = lines.find(target)
    if p > 0:
        return safeReadfloat(lines[p + len(target):p+30])
    else:
        return None

def parsePhiResult(lines):
    """Parse the result of a standard Phi run."""
    target = "PHI (Normal):"
    return findFloatTarget(lines, target)

def parsePhiExtra(lines):
    """Parse the results of calling Phi with the -o option."""
    target1 = "NSS:"
    target2 = "Max Chi^2:"
    target3 = "PHI (Normal):"
    return {'nss': findFloatTarget(lines, target1),
            'chi': findFloatTarget(lines, target2),
            'phi': findFloatTarget(lines, target3)}

def readSequences(filename):
    raw =  SeqIO.parse(filename, "fasta")
    return [r for r in raw]

def countSequences(records):
    n = 0
    for rec in records:
        n += 1
    return n

def findAlignmentLength(records):
    """Return the length of the shortest sequence in `records'."""
    l = 1000000                 # we're probably not going to use this program on very long sequences
    for rec in records:
        l = min(l, len(rec))
    return l

def writeWindowedSequences(filename, records, start=0, end=10):
    with open(filename, "w") as out:
        for rec in records:
            out.write(rec[start:end].format("fasta"))

# Main loop(s)

def doSwipt(G, records, alnlen):
    tmp = tmpfile()
    cmdline = [G.PHI, "-f", tmp]
    end = alnlen - G.WINSIZE
    out = sys.stdout
    if G.OUTFILE:
        out = open(G.OUTFILE, "w")
    # print range(0, end, G.WINSTEP)
    out.write("Start\tEnd\tPval\n")
    try:
        for p in range(0, end, G.WINSTEP):
            q = p+G.WINSIZE
            writeWindowedSequences(tmp, records, p, q)
            result = invokePhi(cmdline)
            pval = parsePhiResult(result)
            if pval == None:
                out.write("{}\t{}\tN/A\n".format(p+1, q))
            elif (not G.PVALTHR) or (pval <= G.PVALTHR):
                out.write("{}\t{}\t{}\n".format(p+1, q, pval))
    finally:
        if G.OUTFILE:
            out.close()
        os.remove(tmp)

def doSwiptExtra(G, records, alnlen):
    tmp = tmpfile()
    cmdline = [G.PHI, "-o", "-f", tmp]
    end = alnlen - G.WINSIZE
    out = sys.stdout
    if G.OUTFILE:
        out = open(G.OUTFILE, "w")
    # print range(0, end, G.WINSTEP)
    out.write("Start\tEnd\tPHI\tNSS\tMax Chi^2\n")
    try:
        for p in range(0, end, G.WINSTEP):
            q = p+G.WINSIZE
            writeWindowedSequences(tmp, records, p, q)
            result = invokePhi(cmdline)
            pvals = parsePhiExtra(result)
            pval = pvals['phi']
            if pval == None:
                out.write("{}\t{}\tN/A\t{}\t{}\n".format(p+1, q, pvals['nss'], pvals['chi']))
            elif (not G.PVALTHR) or (pval <= G.PVALTHR):
                out.write("{}\t{}\t{}\t{}\t{}\n".format(p+1, q, pval, pvals['nss'], pvals['chi']))
    finally:
        if G.OUTFILE:
            out.close()
        os.remove(tmp)
    
# Top-level

def usage(G):
    sys.stderr.write("""Usage: swipt.py [options] infile [outfile]

infile = file containing input sequence alignment in FASTA format
outfile = file results are written to, in tab-delimited format
          (if unspecified, results will be written to standard output).

Options:
  -p PATH      Path to the Phi executable
  -w W         Set window size to W (default: {})
  -s S         Set window step to S (default: {})
  -t T         Set p-value threshold to T
  -phip N      Set Phi permutations to N (Phi -p option)
  -phio        Output extended Phi results (Phi -o option)
""".format(G.WINSIZE, G.WINSTEP))

def parseArguments(G, args):
    next = ""
    for a in args:
        if a == "-phio":
            G.PHIEXTRA = True
        elif next == "-p":
            G.PHI = a
            next = ""
        elif next == "-w":
            G.WINSIZE = int(a)
            next = ""
        elif next == "-s":
            G.WINSTEP = int(a)
            next = ""
        elif next == "-t":
            G.PVALTHR = float(a)
            next = ""
        elif next == "-phip":
            G.PHIPERM = int(a)
            next = ""
        elif a in ["-p", "-w", "-s", "-t", "-phip"]:
            next = a
        elif G.INFILE:
            G.OUTFILE = a
        else:
            G.INFILE = a

def main():
    G = Globals()
    parseArguments(G, sys.argv[1:])
    if G.INFILE:
        records = readSequences(G.INFILE)
        alncnt = countSequences(records)
        alnlen = findAlignmentLength(records)
        sys.stderr.write("swipt.py - Sliding WIndow Phi Test\n\n")
        sys.stderr.write("Read alignment of {} sequences, length={}\n".format(alncnt, alnlen))
        sys.stderr.write("Window size={}, step={}\n".format(G.WINSIZE, G.WINSTEP))
        if G.OUTFILE:
            sys.stderr.write("Writing output to file {}\n".format(G.OUTFILE))
        if G.PHIEXTRA:
            sys.stderr.write("Reporting additional Phi statistics.\n")
            doSwiptExtra(G, records, alnlen)
        else:
            doSwipt(G, records, alnlen)
    else:
        usage(G)

if __name__ == "__main__":
    main()
