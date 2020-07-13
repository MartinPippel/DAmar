# DAmar
Long read QC, assembly and scaffolding pipelines

DAmar, is a hybrid of our earlier [Marvel](https://github.com/schloi/MARVEL), [Dazzler](https://github.com/thegenemyers), and [Daccord](https://gitlab.com/german.tischler/daccord) systems. We aim to create a front-to-end assembler of pacbio or oxford nanopore long-read sequencing data. This pipeline produces a number of QC metrics at various stages as well as incorporating further technologies including Bionano, 10x and HiC data to scaffold the created contigs. The wrapper scripts included here rely on proprietary software (e.g. Bionano Solve) alongside a number of open-source software.
Published code is accurate and runs at time of article submission. The code is currently undergoing refactoring. This will be completed soon.

2020-07-13

# Requirements 
* networkx 2.1+ (Python library)
---
* GTK3 3.x (optional, visualize overlaps with LAexplorer)
* libhdf5 (optional, extract sequences from HDF5 files)


# Compilation
DAmar uses the GNU autoconf/automake tool set. It can be compiled on Linux using:

```
autoreconf -i -f
./configure
make
```

Running autoreconf requires a complete set of tools including autoconf, automake, autoheader, aclocal and libtool.
