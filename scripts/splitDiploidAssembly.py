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

bcount = 0
scount = 0 
hcount = 0

for fn in fIn:
    f=open(fn, 'r')

    if "filtered" in f.name:
        faprefix = f.name.split(".")[1] + "_"
    else:
        faprefix = ""

    for (name, args, seq) in fasta_iter(f):
        ends = args["ends"]
        (e1, e2) = [ int(x) for x in ends ]

        path = args["path"]
        reads = args["reads"]
        sreads = args["sreads"]
        length = args["length"]

        if e1 == e2:
          fout_b.write(">contig_{}{}_{} path={} ends={} length={} reads={} sreads={}".format(faprefix, bcount, len(seq), path[0], ",".join(ends), ",".join(length), ",".join(reads), ",".join(sreads)))
          fout_b.write("\n{}\n".format(wrap_seq(seq, 100)))
          bcount+=1            
          continue

        if e1 != -1 or e2 != -1 and length < 100000:
          fout_s.write(">contig_{}{}_{} path={} ends={} length={} reads={} sreads={}".format(faprefix, scount, len(seq), path[0], ",".join(ends), ",".join(length), ",".join(reads), ",".join(sreads)))
          fout_s.write("\n{}\n".format(wrap_seq(seq, 100)))
          scount+=1            
          continue

        fout_h.write(">contig_{}{}_{} path={} ends={} length={} reads={} sreads={}".format(faprefix, hcount, len(seq), path[0], ",".join(ends), ",".join(length), ",".join(reads), ",".join(sreads)))
        fout_h.write("\n{}\n".format(wrap_seq(seq, 100)))
        hcount+=1        
    f.close()
    
f.close(fout_b)
f.close(fout_s)
f.close(fout_h)
