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
fIn    = sys.argv[1]
nout   = sys.argv[2]

if len(sys.argv) > 3:
    maxLen = sys.argv[3]
else:
    maxLen = 5000

if len(sys.argv) > 4:
    ovh    = sys.argv[4]
else:
    ovh    = 1000
    
if maxLen < ovh + 1500:    
    sys.exit('Maximum read length {} must be at least 1500 bases larger than ovh {}!'.format(maxLen, ovh))
    
# create ouptut file handler 
fout = open(nout, "w")

f=open(fIn, 'r')

c=0
for (name, args, seq) in fasta_iter(f):
    
    numSplits = len(seq)//maxLen
    
    if numSplits:
        splitLastTwoEvenly=0
        if len(seq) % maxLen < ovh:
            splitLastTwoEvenly=1
        beg=0
        for i in range(1,numSplits):
            if i == numSplits and splitLastTwoEvenly:
                end = beg + (len(seq) - beg)//2
            else:
                end=beg+maxLen
            fout.write(">splitRead/{}/{}_{} oldName={} beg={}".format(c,beg,end,name))
            fout.write("\n{}\n".format(wrap_seq(seq[beg:end], 100)))
            beg=end-ovh
            c=c+1
        fout.write(">splitRead/{}/{}_{} oldName={} beg={}".format(c,beg,len(seq),name))
        fout.write("\n{}\n".format(wrap_seq(seq[beg:len(seq)], 100)))
        c=c+1        
    else:
        fout.write(">splitRead/{}/{}_{} oldName={} beg={}".format(c,0,len(seq),name))
        fout.write("\n{}\n".format(wrap_seq(seq, 100)))
        c=c+1            
f.close()
fout.close()    

