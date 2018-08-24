#!/usr/bin/env python


import os
import itertools
import sys
import string
import operator

from __future__ import print_function
import sys

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

ng_size = float(sys.argv[1]) 
if len(sys.argv) > 2:
    filein = sys.argv[2]
else:
    filein = None

if filein == None or filein == "-":
    fIn = sys.stdin
else:
    fIn = open(filein)

arrSizes = []
for (name, args, seq) in fasta_iter(fIn):
    arrSizes.append( (name, args, len(seq), seq) )

arrSizes.sort( key = operator.itemgetter(2), reverse = True )

n_size  = sum( x[2] for x in arrSizes )
n_total = 0
n_count = 0

while n_total < n_size / 2:
    n_total += arrSizes[n_count][2]
    n_count += 1

ng_total = 0
ng_count = 0

if n_size < ng_size:
    eprint("Assembly size {} is smaller than estimated genome size {}".format(n_size, ng_size))
if n_size < ng_size/2:
    eprint("Assembly size {} is smaller than given NG50 {}. Reset to N50.".format(n_size, ng_size/2))
    ng_size = n_size
    ng_total= n_total
    ng_count= n_count
else:
    while gtotal < ng_size / 2:
        ng_total += arrSizes[gcount][2]
        ng_count += 1

i = 1
nSmallContigs=0
n1MContigs=0
n5MContigs=0
n10MContigs=0
n20MContigs=0
nlt1MContigs=0

for i in range( len(arrSizes) ):
    if arrSizes[i][2] < 200000:
       nSmallContigs+=1
    elif arrSizes[i][2] >= 20000000:
       n20MContigs+=1
    elif arrSizes[i][2] >= 10000000:
       n10MContigs+=1
    elif arrSizes[i][2] >= 5000000:
       n5MContigs+=1
    elif arrSizes[i][2] >= 1000000:
       n1MContigs+=1
    else:
       nlt1MContigs+=1

cumSeq=0
firstNg10=1
firstNg20=1
firstNg30=1
firstNg40=1
firstNg50=1
firstNg60=1
firstNg70=1
firstNg80=1
firstNg90=1
firstNg100=1

for i in range( len(arrSizes) ):
#    reads = arrSizes[i][1]["reads"][::3]
#    sreads = arrSizes[i][1]["sreads"]
#    source = ",".join( "{}-{}".format(*x) for x in zip(reads, sreads) )

    seq = arrSizes[i][3]
    cumSeq += len(seq)
    #print("{:2} {:8} {:10}".format(i, arrSizes[i][2], cumSeq))
    sys.stdout.write("{:2} {:8} {:10}".format(i, arrSizes[i][2], cumSeq))
    if firstNg10 == 1 and cumSeq >= 0.1*ng_size:
      sys.stdout.write(" NG10")
      firstNg10=0
    if firstNg20 == 1 and cumSeq >= 0.2*ng_size:
      sys.stdout.write(" NG20")
      firstNg20=0
    if firstNg30 == 1 and cumSeq >= 0.3*ng_size:
      sys.stdout.write(" NG30")
      firstNg30=0
    if firstNg40 == 1 and cumSeq >= 0.4*ng_size:
      sys.stdout.write(" NG40")
      firstNg40=0
    if firstNg50 == 1 and cumSeq >= 0.5*ng_size:
      sys.stdout.write(" NG50")
      firstNg50=0
    if firstNg60 == 1 and cumSeq >= 0.6*ng_size:
      sys.stdout.write(" NG60")
      firstNg60=0
    if firstNg70 == 1 and cumSeq >= 0.7*ng_size:
      sys.stdout.write(" NG70")
      firstNg70=0
    if firstNg80 == 1 and cumSeq >= 0.8*ng_size:
      sys.stdout.write(" NG80")
      firstNg80=0
    if firstNg90 == 1 and cumSeq >= 0.9*ng_size:
      sys.stdout.write(" NG90")
      firstNg90=0
    if firstNg100 == 1 and cumSeq >= 1.0*ng_size:
      sys.stdout.write(" NG100")
      firstNg100=0    
    sys.stdout.write("\n")
 
print("{} in {} ... ng50 {} ... asm size {}".format(ng_total, ng_count, arrSizes[ng_gcount-1][2], ng_size))
print("{} in {} ... n50 {}  ... asm size {}".format(n_total, n_count, arrSizes[n_count-1][2], n_size))
print("longest_3 {} {} {}".format(arrSizes[0][2], arrSizes[1][2], arrSizes[2][2]))
print("contigs lt200k {} lt1M {} ge1M {} ge5M {} ge10M {} ge20M {}".format(nSmallContigs,nlt1MContigs,n1MContigs,n5MContigs,n10MContigs,n20MContigs))
