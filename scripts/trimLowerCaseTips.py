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
        header = next(header)[1:].strip()
        seq = "".join(s.strip() for s in next(faiter))

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
fIn= sys.argv[1]
nout= sys.argv[2]
nout_garbage= sys.argv[3]

# create ouptut file handler 
fout = open(nout, "w")
fout_garbage = open(nout_garbage, "w")

f=open(fIn, 'r')

for (name, args, seq) in fasta_iter(f):
    
    slen   =len(seq)
    trimBeg=0
    trimEnd=0
    
    for i in seq:
        if i.islower():
            trimBeg+=1
        else:
            break
        
    for i in reversed(seq):
        if i.islower():
            trimEnd+=1
        else:
            break
        
    if trimEnd > 0 and trimEnd < slen:
        trimEnd=-trimEnd
        
    if trimEnd == 0:
        trimEnd = slen
    
    if trimBeg == slen:
        fout_garbage.write(">{} trimBeg={} trimEnd={}\n".format(name, trimBeg, trimEnd))
        fout_garbage.write("\n{}\n".format(wrap_seq(seq, 100)))
    else:
        fout.write(">{} trimBeg={} trimEnd={}".format(name, trimBeg, slen - abs(trimEnd)))
        fout.write("\n{}\n".format(wrap_seq(seq[trimBeg:trimEnd], 100)))
        
f.close()
fout.close()
fout_garbage.close()    

