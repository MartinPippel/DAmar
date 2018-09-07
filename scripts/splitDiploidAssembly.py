#!/usr/bin/env python
from __future__ import print_function
import glob
import os
import itertools
import sys 
import string
import operator

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

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
fout_p = open("{}/{}.p.fasta".format(outDir,prefix), "w")
fout_a = open("{}/{}.a.fasta".format(outDir,prefix), "w")

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
          fout_a.write(">{} path={} ends={} length={} reads={} sreads={}".format(name, path[0], ",".join(ends), length, ",".join(reads), ",".join(sreads)))
          fout_a.write("\n{}\n".format(wrap_seq(seq, 100)))
          if (length > 500000):
              eprint("WARNING Found very long alternative Contig (bubble): {} lenght: {}.".format(name, length))
          continue

        if (e1 != -1 or e2 != -1) and length < 100000:
          fout_a.write(">{} path={} ends={} length={} reads={} sreads={}".format(name, path[0], ",".join(ends), length, ",".join(reads), ",".join(sreads)))
          fout_a.write("\n{}\n".format(wrap_seq(seq, 100)))
          continue

        fout_p.write(">{} path={} ends={} length={} reads={} sreads={}".format(name, path[0], ",".join(ends), length, ",".join(reads), ",".join(sreads)))
        fout_p.write("\n{}\n".format(wrap_seq(seq, 100)))
    f.close()
    
fout_p.close()
fout_a.close()
