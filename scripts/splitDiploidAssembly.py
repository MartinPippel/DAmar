#!/usr/bin/env python

import glob
import os
import itertools
import sys 
import string
import operator

def fasta_iter(fh):
    faiter = (x[1] for x in itertools.groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        header = header.next()[1:].strip()
        seq = "".join(s.strip() for s in faiter.next())

        fields = header.split()
        args = {}
        for field in fields[1:]:
            (name, dummy, values) = field.partition("=")
            args[name] = values.split(",")

        yield ( fields[0], args, seq )

def wrap_seq(strSeq, strLen):
    arrLines = []

    for i in range(0, len(strSeq), strLen):
        arrLines.append( strSeq[i : i+strLen] )

    return "\n".join(arrLines)

# outfile prefix
prefix = sys.argv[1]
nsize = float(sys.argv[2]) 
outDir = sys.argv[3]

if len(sys.argv) > 4:
    filein = sys.argv[4]
else:
    filein = None

if filein == None or filein == "-":
    fIn = [ sys.stdin ]
else:
#    if os.path.isdir(filein):
#        fIn = [ open(x) for x in glob.glob( os.path.join(filein, "*.fasta") ) ]
#    else:
#        fIn = [ open(filein) ]
    if os.path.isdir(filein):
        fIn = glob.glob( os.path.join(filein, "*.fasta") ) 
    else:
        fIn = [ filein ]

# create ouptut file handler 
fout_h = open("{}/{}.haploid.fasta".format(outDir,prefix), "w")
fout_b = open("{}/{}.bubbles.fasta".format(outDir,prefix), "w")
fout_s = open("{}/{}.spurs.fasta".format(outDir,prefix), "w")

for fn in fIn:
    f=open(fn, 'r')

    for (name, args, seq) in fasta_iter(f):
        ends = args["ends"]
        (e1, e2) = [ int(x) for x in ends ]

        path = args["path"]
        reads = args["reads"]
        sreads = args["sreads"]
        length = int(args["length"][0])

        if e1 == e2 and e1 != -1:
          fout_b.write(">{} path={} ends={} length={} reads={} sreads={}".format(name, path[0], ",".join(ends), length, ",".join(reads), ",".join(sreads)))
          fout_b.write("\n{}\n".format(wrap_seq(seq, 100)))
          continue

        if (e1 != -1 or e2 != -1) and length < 100000:
          fout_s.write(">{} path={} ends={} length={} reads={} sreads={}".format(name, path[0], ",".join(ends), length, ",".join(reads), ",".join(sreads)))
          fout_s.write("\n{}\n".format(wrap_seq(seq, 100)))
          continue

        fout_h.write(">{} path={} ends={} length={} reads={} sreads={}".format(name, path[0], ",".join(ends), length, ",".join(reads), ",".join(sreads)))
        fout_h.write("\n{}\n".format(wrap_seq(seq, 100)))
    f.close()
    
fout_b.close()
fout_s.close()
fout_h.close()
