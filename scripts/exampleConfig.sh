####### MARVEL pipeline variable 

### MARVEL PATH
MARVEL_SOURCE_PATH="/projects/dazzler/pippel/prog/MARVEL/DAMARVEL"
MARVEL_PATH="/projects/dazzler/pippel/prog/MARVEL/DAMARVEL-build"

### REPCOMP PATH
REPCOMP_SOURCE_PATH="/projects/dazzler/pippel/prog/repcomp"
REPCOMP_PATH="/projects/dazzler/pippel/prog/repcomp-build"

### DACCORD PATH - used progs: fastaidrename, forcealign
DACCORD_SOURCE_PATH="/projects/dazzler/pippel/prog/daccord/"
DACCORD_PATH="/projects/dazzler/pippel/prog/daccord-build/"

### DZZLER PATH
DAZZLER_SOURCE_PATH="/projects/dazzler/pippel/prog/dazzler/"
DAZZLER_PATH="/projects/dazzler/pippel/prog/dazzler/"

### slurm scripts path
SUBMIT_SCRIPTS_PATH="${MARVEL_PATH}/scripts"

############################## tools for pacbio arrow correction 
### samtools
SAMTOOLS_PATH="/projects/dazzler/pippel/prog/samtools/samtools-1.8"
### bamtools 
BAMTOOLS_PATH="/projects/dazzler/pippel/prog/bamtools/build"
### arrow path
PACBIO_ARROW_TOOLS="/projects/dazzler/pippel/prog/pacbio/FALCON_ALL_PREBUILD"

### ENVIRONMENT VARIABLES 
export PYTHONUSERBASE=${PACBIO_ARROW_TOOLS}
export LD_LIBRARY_PATH=${PACBIO_ARROW_TOOLS}/lib:${LD_LIBRARY_PATH}
export PATH=${PACBIO_ARROW_TOOLS}/bin:${PACBIO_ARROW_TOOLS}:${PATH}
export PATH=${MARVEL_PATH}/bin:${MARVEL_PATH}/scripts:$PATH
export PYTHONPATH=${MARVEL_PATH}/lib.python:$PYTHONPATH
export PATH=${SAMTOOLS_PATH}/bin:$PATH
export PATH=${BAMTOOLS_PATH}/bin:$PATH

## general information
PROJECT_ID=iHylVes1
GSIZE=600M
SLURM_PARTITION=batch

## general settings raw read data bases  
RAW_DB=LAB1608_HYLES_VESPERTILIO_MARVEL
RAW_DAZZ_DB=LAB1608_HYLES_VESPERTILIO_DAZZLER
RAW_COV=40

## general settings patched read databases 
FIX_DB=LAB1608_HYLES_VESPERTILIO_MARVEL_FIX
FIX_DAZZ_DB=LAB1608_HYLES_VESPERTILIO_DAZZLER_FIX
FIX_COV=40

## corrected DB 
COR_DB=LAB1608_HYLES_VESPERTILIO_MARVEL_COR

## corrected contig DB
CONT_DB=LAB1608_HYLES_VESPERTILIO_MARVEL_CONTIG
CONT_DAZZ_DB=LAB1608_HYLES_VESPERTILIO_DAZZLER_CONTIG
 
################# define marvel phases and their steps that should be done 

DB_PATH=/projects/dazzlerAssembly/LAB1608.HYLES_VESPERTILIO/data/pacbio
PATCHING_DIR="patching"
ASSMEBLY_DIR="assembly"

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> marvel phase 1 - repeat masking <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# type-0 steps [1-13]: 1-DBdust, 2-Catrack, 3-datander, 4-TANmask, 5-Catrack, 6-daligner, 7-LAmerge, 8-LArepeat, 9-TKmerge, 10-daligner, 11-LAmerge, 12-LArepeat, 13-TKmerge]
RAW_REPMASK_TYPE=0

RAW_REPMASK_SUBMIT_SCRIPTS_FROM=1
RAW_REPMASK_SUBMIT_SCRIPTS_TO=9

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> marvel phase 2 - read patching <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 
#type-0 steps   [1-11]:   1-daligner, 2-LAmerge, 3-LArepeat, 4-TKmerge, 5-TKcombine, 6-TKhomogenize, 7-TKcombine, 8-LAfilter, 9-LAq, 10-TKmerge, 11-LAfix                 #old pipeline
#type-1 steps   [ 1-11 :  1-daligner, 2-LAmerge, 3-LArepeat, 4-TKmerge, 5-TKcombine, 6-TKhomogenize, 7-TKcombine, 8-LAfilter, 9-LAq, 10-TKmerge, 11-LAfix                # experimental pipeline     
#                12-22 : 12-LAseparate, 13-repcomp, 14-LAmerge, 15-LArepeat, 16-TKmerge, 17-TKcombine, 18-TKhomogenize, 19-TKcombine, 21-LAq, 21-TKmerge, 22-LAfix
#                23-33]: 23-LAseparate, 24-forcealign, 25-LAmerge, 26-LArepeat, 27-TKmerge, 28-TKcombine, 29-TKhomogenize, 30-TKcombine, 31-LAq, 32-TKmerge, 33-LAfix
RAW_PATCH_TYPE=1

RAW_PATCH_SUBMIT_SCRIPTS_FROM=1
RAW_PATCH_SUBMIT_SCRIPTS_TO=33

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> marvel phase 3 - repeat masking <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# type_0 steps [1-14]: 1-createFIX_DB, 2-DBdust, 3-Catrack, 4-datander, 5-TANmask, 6-Catrack, 7-daligner, 8-LAmerge, 9-LArepeat, 10-TKmerge, 11-daligner, 12-LAmerge, 13-LArepeat, 14-TKmerge]
FIX_REPMASK_TYPE=0

FIX_REPMASK_SUBMIT_SCRIPTS_FROM=1
FIX_REPMASK_SUBMIT_SCRIPTS_TO=10

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> marvel phase 4 - scrubbing <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
FIX_SCRUB_TYPE=1

#type-0 steps [1-13]: 1-daligner, 2-LAmerge, 3-LArepeat, 4-TKmerge, 5-TKcombine, 6-TKhomogenize, 7-TKcombine, 8-LAstitch, 9-LAq, 10-TKmerge, 11-LAgap, 12-LAq, 13-TKmerge          ## old pipeline
#type-1 steps [ 1-13 : 1-daligner, 2-LAmerge, 3-LArepeat, 4-TKmerge, 5-TKcombine, 6-TKhomogenize, 7-TKcombine, 8-LAstitch, 9-LAq, 10-TKmerge, 11-LAgap, 12-LAq, 13-TKmerge          ## experimental pipeline
#              14-27 : 14-LAseparate, 15-repcomp, 16-LAmerge, 17-LArepeat, 18-TKmerge, 19-TKcombine, 20-TKhomogenize, 21-TKcombine, 22-LAstitch, 23-LAq, 24-TKmerge, 25-LAgap, 26-LAq, 27-TKmerge
#              28-41]: 28-LAseparate, 29-forcealign, 30-LAmerge, 31-LArepeat, 32-TKmerge, 33-TKcombine, 34-TKhomogenize, 35-TKcombine, 36-LAstitch, 37-LAq, 38-TKmerge, 39-LAgap, 40-LAq, 41-TKmerge
FIX_SCRUB_SUBMIT_SCRIPTS_FROM=1
FIX_SCRUB_SUBMIT_SCRIPTS_TO=41

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> marvel phase 5 - overlap filtering <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

FIX_FILT_TYPE=0
#type-1 steps: 1-createSubdir, 2-LAfilter, 3-LAmerge

FIX_FILT_SUBMIT_SCRIPTS_FROM=1
FIX_FILT_SUBMIT_SCRIPTS_TO=3

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> marvel phase 6 - touring, stats, layout <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

FIX_TOUR_TYPE=0
#type-0 steps: 1-OGbuild, 2-OGtour, 3-tour2fasta, 4-OGlayout, 5-statistics
FIX_TOUR_SUBMIT_SCRIPTS_FROM=1
FIX_TOUR_SUBMIT_SCRIPTS_TO=5

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> marvel phase 7 - contig correction  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

FIX_CORR_TYPE=0
#type-0 steps: 1-paths2rids, 2-LAcorrect, 3-prepDB, 4-tour2fasta, 5-statistics
FIX_CORR_SUBMIT_SCRIPTS_FROM=1
FIX_CORR_SUBMIT_SCRIPTS_TO=5

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> marvel phase 8 - contig analyze stats  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

COR_CONTIG_TYPE=0
#type-0 steps: 1-createCorrectedContigDB, 2-DBdust, 3-Catrack, 4-datander, 5-TANmask, 6-Catrack, 7-daligner, 8-LAfilter, 9-LAseparate, 10-forcealign, 11-LAmerge, 12-LArepeat, 13-TKmerge, 14-TKcombine, 15-LAfilter, 16-LAmerge, 17-CTanalyze
COR_CONTIG_SUBMIT_SCRIPTS_FROM=0
COR_CONTIG_SUBMIT_SCRIPTS_TO=17

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> marvel phase 9 - pacbio arrow correction  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

PB_ARROW_TYPE=0
#type-0 steps: 1-prepInFasta, 2-pbalign, 3-bamsplit, 4-bamseparate, 5-bamMerge, 6-arrow, 7-statistics
PB_ARROW_SUBMIT_SCRIPTS_FROM=1
PB_ARROW_SUBMIT_SCRIPTS_TO=7


# ----------------------------------------------------------------- RAW REPEAT MASKING OPTIONS - always on RAW_DB ---------------------------------------------------------------------------

# number of daligner block comparison for repeat masking if there are more then one number then mulyiple round of repeat masking are applied
RAW_REPMASK_BLOCKCMP=(2)
# bias in nucleotide distribution - set to 1 for AT >= 70% || AT <= 30% 
RAW_REPMASK_DBDUST_BIAS=0
# datander
RAW_REPMASK_DATANDER_THREADS=8
RAW_REPMASK_DATANDER_MINLEN=500
RAW_REPMASK_DATANDER_FOLDER="tan"
# TANmask
RAW_REPMASK_TANMASK_TRACK="rawtan"
RAW_REPMASK_TANMASK_MINLEN=500
RAW_REPMASK_TANMASK_VERBOSE=1
# catrack 
RAW_REPMASK_CATRACK_DELETE=1
RAW_REPMASK_CATRACK_VERBOSE=1
RAW_REPMASK_CATRACK_OVERWRITE=1
# daligner
RAW_REPMASK_DALIGNER_RUNID=0
RAW_REPMASK_DALIGNER_VERBOSE=1
RAW_REPMASK_DALIGNER_BIAS=${REPMASK_DBDUST_BIAS}
RAW_REPMASK_DALIGNER_IDENTITY_OVLS=0
RAW_REPMASK_DALIGNER_KMER=14
RAW_REPMASK_DALIGNER_OLEN=700
RAW_REPMASK_DALIGNER_ERR=0.7
RAW_REPMASK_DALIGNER_MASK="dust ${RAW_REPMASK_TANMASK_TRACK}"
# LArepeat
RAW_REPMASK_LAREPEAT_COV=(20)
RAW_REPMASK_LAREPEAT_LOW=1
RAW_REPMASK_LAREPEAT_HGH=1
RAW_REPMASK_LAREPEAT_REPEATTRACK=rawRepmask
# TKmerge 
RAW_REPMASK_TKMERGE_DELETE=1	# i.e. delete intermediate block based tracks

# ----------------------------------------------------------------- RAW PATCHING OPTIONS  --------------------------------------------------------------------------------------------------- 

# daligner
RAW_FIX_DALIGNER_RUNID=$((${RAW_REPMASK_DALIGNER_RUNID}+1))
RAW_FIX_DALIGNER_VERBOSE=1
RAW_FIX_DALIGNER_DAL=22 
RAW_FIX_DALIGNER_BIAS=${REPMASK_DBDUST_BIAS}
RAW_FIX_DALIGNER_IDENTITY_OVLS=1
x=0
RAW_FIX_DALIGNER_MASK="dust ${RAW_REPMASK_TANMASK_TRACK}"
while [[ $x -lt ${#RAW_REPMASK_BLOCKCMP[*]} ]]
do 
	RAW_FIX_DALIGNER_MASK="${RAW_FIX_DALIGNER_MASK} ${RAW_REPMASK_LAREPEAT_REPEATTRACK}_B${RAW_REPMASK_BLOCKCMP[${x}]}C${RAW_REPMASK_LAREPEAT_COV[${x}]}"
	x=$(($x+1))
done
RAW_FIX_DALIGNER_NUMACTL=1
#  LAmerge
RAW_FIX_LAMERGE_NFILES=64
# LArepeat
RAW_FIX_LAREPEAT_LEAVE_COV=1.7
RAW_FIX_LAREPEAT_ENTER_COV=2.0
RAW_FIX_LAREPEAT_COV=${RAW_COV}
# TKmerge
RAW_FIX_TKMERGE_DELETE=1
# TKcombine
RAW_FIX_TKCOMBINE_DELETE=1
RAW_FIX_TKCOMBINE_VERBOSE=1 
# LAseparate
RAW_FIX_LASEPARATE_USEREPEAT=1 
# repcomp
RAW_FIX_REPCOMP_INBLOCKSIZE=64M
RAW_FIX_REPCOMP_MEM=64g
RAW_FIX_REPCOMP_KMER=14
RAW_FIX_REPCOMP_THREADS=4
RAW_FIX_REPCOMP_RUNID=$((${RAW_FIX_DALIGNER_RUNID}+1))
RAW_FIX_REPCOMP_NUMACTL=1
RAW_FIX_REPCOMP_CORRELATION=0.7
RAW_FIX_REPCOMP_MASK="dust"
RAW_FIX_REPCOMP_OLEN=700
# forcealign 
RAW_FIX_FORCEALIGN_THREADS=1
RAW_FIX_FORCEALIGN_RUNID=$((${RAW_FIX_REPCOMP_RUNID}+1))
RAW_FIX_FORCEALIGN_PARTIAL=1
RAW_FIX_FORCEALIGN_NUMACTL=1
RAW_FIX_FORCEALIGN_MAXDIST=250
RAW_FIX_FORCEALIGN_BORDER=25
RAW_FIX_FORCEALIGN_CORRELATION=0.65
# LAq
RAW_FIX_LAQ_MINSEG=1
RAW_FIX_LAQ_QTRIMCUTOFF=35 # thats always related to a 100 bp trace spacing and will adjusted appropriately if it differs 
# LAfix
#RAW_FIX_LAFIX_QUAL=30    # default is 28, only increase if data has very low quality 
RAW_FIX_LAFIX_GAP=-1
RAW_FIX_LAFIX_MLEN=4000
RAW_FIX_LAFIX_LOW_COVERAGE=0 
RAW_FIX_LAFIX_TRIM=1
#RAW_FIX_LAFIX_QTRACK=q0_d35_s1 ## will be set automatically according to LAq settings
RAW_FIX_LAFIX_USEREPEAT=1
RAW_FIX_LAFIX_FILESUFFIX=".patchType${RAW_PATCH_TYPE}"
RAW_FIX_LAFIX_CONVERTRACKS=""
while [[ $x -lt ${#RAW_REPMASK_BLOCKCMP[*]} ]]
do 
	RAW_FIX_LAFIX_CONVERTRACKS="${RAW_FIX_LAFIX_CONVERTRACKS} ${RAW_REPMASK_LAREPEAT_REPEATTRACK}_B${RAW_REPMASK_BLOCKCMP[${x}]}C${RAW_REPMASK_LAREPEAT_COV[${x}]}"
	x=$(($x+1))
done 
RAW_FIX_LAFIX_DISCARDCHIMERS=1
RAW_FIX_LAFIX_MINCHIMERBORDERCOV=7
RAW_FIX_LAFIX_MAXCHIMERLEN=8000

# ----------------------------------------------------------------- FIX REPEAT MASKING OPTIONS - always on RAW_DB ---------------------------------------------------------------------------

FIX_REPMASK_USELAFIX_PATH=patchedReads_dalign
# number of daligner block comparison for repeat masking if there are more then one number then mulyiple round of repeat masking are applied
FIX_REPMASK_BLOCKCMP=(2)
# bias in nucleotide distribution - set to 1 for AT >= 70% || AT <= 30% 
FIX_REPMASK_DBDUST_BIAS=0
# datander
FIX_REPMASK_DATANDER_THREADS=8
FIX_REPMASK_DATANDER_MINLEN=500
FIX_REPMASK_DATANDER_FOLDER="tan"
# TANmask
FIX_REPMASK_TANMASK_TRACK="tan"
FIX_REPMASK_TANMASK_MINLEN=500
FIX_REPMASK_TANMASK_VERBOSE=1
# catrack 
FIX_REPMASK_CATRACK_DELETE=1
FIX_REPMASK_CATRACK_VERBOSE=1
FIX_REPMASK_CATRACK_OVERWRITE=1
# daligner
FIX_REPMASK_DALIGNER_RUNID=0
FIX_REPMASK_DALIGNER_VERBOSE=1
FIX_REPMASK_DALIGNER_BIAS=${REPMASK_DBDUST_BIAS}
FIX_REPMASK_DALIGNER_IDENTITY_OVLS=1
FIX_REPMASK_DALIGNER_KMER=14
FIX_REPMASK_DALIGNER_OLEN=700
FIX_REPMASK_DALIGNER_ERR=0.7
FIX_REPMASK_DALIGNER_MASK="dust ${FIX_REPMASK_TANMASK_TRACK}"
# LArepeat
FIX_REPMASK_LAREPEAT_COV=(20)
FIX_REPMASK_LAREPEAT_LOW=1
FIX_REPMASK_LAREPEAT_HGH=1
FIX_REPMASK_LAREPEAT_REPEATTRACK=repmask
# TKmerge 
FIX_REPMASK_TKMERGE_DELETE=1	# i.e. delete intermediate block based tracks

# ----------------------------------------------------------------- SCRUBBING OPTIONS FOR PATCHED DB ----------------------------------------------------------------------------------------

# daligner
FIX_SCRUB_DALIGNER_RUNID=$((${FIX_REPMASK_DALIGNER_RUNID}+1))
FIX_SCRUB_DALIGNER_DAL=12 # number of block comparison per daligner job, to keep somehoew below 1000 jobs 
FIX_SCRUB_DALIGNER_BIAS=${FIX_REPMASK_DBDUST_BIAS}
x=0
FIX_SCRUB_DALIGNER_MASK="dust ${FIX_REPMASK_TANMASK_TRACK}"
while [[ $x -lt ${#FIX_REPMASK_BLOCKCMP[*]} ]]
do 
	FIX_SCRUB_DALIGNER_MASK="${FIX_SCRUB_DALIGNER_MASK} ${FIX_REPMASK_LAREPEAT_REPEATTRACK}_B${FIX_REPMASK_BLOCKCMP[${x}]}C${FIX_REPMASK_LAREPEAT_COV[${x}]}"
	x=$(($x+1))
done 
FIX_SCRUB_DALIGNER_NUMACTL=1
# repcomp
FIX_SCRUB_REPCOMP_INBLOCKSIZE=64M
FIX_SCRUB_REPCOMP_MEM=64g
FIX_SCRUB_REPCOMP_KMER=14
FIX_SCRUB_REPCOMP_THREADS=4
FIX_SCRUB_REPCOMP_RUNID=$((${FIX_SCRUB_DALIGNER_RUNID}+1))
FIX_SCRUB_REPCOMP_NUMACTL=1
FIX_SCRUB_REPCOMP_CORRELATION=0.7
FIX_SCRUB_REPCOMP_MASK="${FIX_REPMASK_TANMASK_TRACK}_dust"
FIX_SCRUB_REPCOMP_OLEN=700
# forcealign 
FIX_SCRUB_FORCEALIGN_THREADS=1
FIX_SCRUB_FORCEALIGN_RUNID=$((${FIX_SCRUB_REPCOMP_RUNID}+1))
FIX_SCRUB_FORCEALIGN_PARTIAL=1
FIX_SCRUB_FORCEALIGN_NUMACTL=1
FIX_SCRUB_FORCEALIGN_MAXDIST=250
FIX_SCRUB_FORCEALIGN_BORDER=25
FIX_SCRUB_FORCEALIGN_CORRELATION=0.65
# LAseparate
FIX_SCRUB_LASEPARATE_USEREPEATIDX=3
# LAmerge
FIX_SCRUB_LAMERGE_NFILES=32
# LArepeat
FIX_SCRUB_LAREPEAT_LEAVE_COV=(1.5 1.7 1.5 1.7)
FIX_SCRUB_LAREPEAT_ENTER_COV=(2.0 2.0 2.0 2.0)
FIX_SCRUB_LAREPEAT_COV=(-1 -1 ${FIX_COV} ${FIX_COV})
# TKmerge 
FIX_SCRUB_TKMERGE_DELETE=1
# TKcombine 
FIX_SCRUB_TKCOMBINE_DELETE=1
FIX_SCRUB_TKCOMBINE_VERBOSE=1 
# stitch
FIX_SCRUB_LASTITCH_FUZZ=100
FIX_SCRUB_LASTITCH_ANCHOR=1000
FIX_SCRUB_LASTITCH_PRELOAD=1
FIX_SCRUB_LASTITCH_PURGE=0
FIX_SCRUB_LASTITCH_LOWCOMPLEXITY=${FIX_REPMASK_TANMASK_TRACK}_dust
FIX_SCRUB_LASTITCH_VERBOSE=0
FIX_SCRUB_LASTITCH_REPEATIDX=0
# LAq
FIX_SCRUB_LAQ_MINSEG=1
FIX_SCRUB_LAQ_QTRIMCUTOFF=27 # thats always related to a 100 bp trace spacing and will adjusted appropriately if it differs 
# LAgap
#FIX_SCRUB_LAGAP_STITCH=100
FIX_SCRUB_LAGAP_PRELOAD=1
FIX_SCRUB_LAGAP_PURGE=0
FIX_SCRUB_LAGAP_TRIM=1
FIX_SCRUB_LAGAP_DISCARD_CHIMERS=0  # repeat defined by index and FIX_SCRUB_LAREPEAT_*

# ----------------------------------------------------------------- FILTER OPTIONS ----------------------------------------------------------------------------------------------------------

# general setting
FIX_FILT_SCRUB_TYPE=1
FIX_FILT_OUTDIR=m1
# LAfilter
FIX_FILT_LAFILTER_VERBOSE=0	  # be carful if enabnled creates huge log files 
#FIX_FILT_LAFILTER_DIF=40
FIX_FILT_LAFILTER_NREP=100
FIX_FILT_LAFILTER_OLEN=4000
FIX_FILT_LAFILTER_PURGE=1
#FIX_FILT_LAFILTER_STITCH=${SCRUB_LAGAP_STITCH}
FIX_FILT_LAFILTER_REPEAT_IDX=0
FIX_FILT_LAFILTER_TRIM=1
FIX_FILT_LAFILTER_UBAS=0
FIX_FILT_LAFILTER_PRELOAD=1
#FIX_FILT_LAFILTER_MERGEREPEATS=100
#FIX_FILT_LAFILTER_MERGEREPEATTIPS=1000
#FIX_FILT_LAFILTER_MINTIPCOV=3
#FIX_FILT_LAFILTER_RMSYMROUNDS=3
#FIX_FILT_LAFILTER_EXCLUDEREADS=blacklist.txt 
#FIX_FILT_LAFILTER_RESOLVE_REPEATS=1         # define greediness  1: m, 2: mm, 3: mmm
#FIX_FILT_LAFILTER_RESOLVE_REPEATS_AGG=1     # even more greedy: exchange m with M
# LAmerge
FIX_FILT_LAMERGE_NFILES=64

# ----------------------------------------------------------------- TOURING OPTIONS ---------------------------------------------------------------------------------------------------------

# OGbuild
FIX_TOUR_OGBUILD_CONT=0
FIX_TOUR_OGBUILD_SPLIT=1
FIX_TOUR_OGBUILD_TRIM=1
# OGtour
FIX_TOUR_OGTOUR_CIRCULAR=1
FIX_TOUR_OGTOUR_LOOKAHAED=6
FIX_TOUR_OGTOUR_DROPINV=1  
# tour2fasta
FIX_TOUR_2FASTA_TRIM=1
FIX_TOUR_2FASTA_SPLIT=0 ### split at every junction
# OGlayout
FIX_TOUR_OGLAYOUT_VERBOSE=1
FIX_TOUR_OGLAYOUT_DIST=200
FIX_TOUR_OGLAYOUT_RMREVERSEEDGE=1

# ----------------------------------------------------------------- CORRECTION OPTIONS ------------------------------------------------------------------------------------------------------
COR_DIR=correction
#path2readids
FIX_CORR_PATHS2RIDS_FILE=${COR_DB%.db}.tour.rids
#LAcorrect
FIX_CORR_LACORRECT_VERBOSE=1
FIX_CORR_LACORRECT_THREAD=1
#tour2fasta
FIX_CORR_2FASTA_SPLIT=0


# ----------------------------------------------------------------- CONTIG ANALYZE OPTIONS --------------------------------------------------------------------------------------------------

ANALYZE_DIR=analyze

### DBdust
COR_CONTIG_DBDUST_BIAS=0
### catrack
COR_CONTIG_CATRACK_VERBOSE=1
COR_CONTIG_CATRACK_DELETE=1
COR_CONTIG_CATRACK_OVERWRITE=1
### datander
COR_CONTIG_DATANDER_THREADS=4
COR_CONTIG_DATANDER_MINLEN=500
COR_CONTIG_DATANDER_FOLDER="tan"
### TANmask 
COR_CONTIG_TANMASK_VERBOSE=1
COR_CONTIG_TANMASK_MINLEN=500
COR_CONTIG_TANMASK_TRACK="tan"
### daligner
COR_CONTIG_DALIGNER_IDENTITY_OVLS=1
COR_CONTIG_DALIGNER_KMER=16
COR_CONTIG_DALIGNER_ERR=0.7
COR_CONTIG_DALIGNER_BIAS=0
COR_CONTIG_DALIGNER_RUNID=1
COR_CONTIG_DALIGNER_OLEN=0
COR_CONTIG_DALIGNER_MEM=64
COR_CONTIG_DALIGNER_MASK="dust ${COR_CONTIG_DATANDER_FOLDER}"
COR_CONTIG_DALIGNER_TRACESPACE=0
### LAmerge
COR_CONTIG_LAMERGE_NFILES=255
### LAseparate
COR_CONTIG_LASEPARATE_OLEN=0
COR_CONTIG_LASEPARATE_RLEN=0
COR_CONTIG_LASEPARATE_REPEAT="${COR_CONTIG_TANMASK_TRACK}"
### forcealign 
COR_CONTIG_FORCEALIGN_THREADS=1
COR_CONTIG_FORCEALIGN_RUNID=$((${COR_CONTIG_DALIGNER_RUNID}+1))
COR_CONTIG_FORCEALIGN_PARTIAL=1
COR_CONTIG_FORCEALIGN_NUMACTL=1
COR_CONTIG_FORCEALIGN_MAXDIST=250
COR_CONTIG_FORCEALIGN_BORDER=25
COR_CONTIG_FORCEALIGN_CORRELATION=0.65
### LArepeat 
COR_CONTIG_LAREPEAT_COV=(2 3 4 5 50 100 1000)
COR_CONTIG_LAREPEAT_ENTER_COV=(1 1 1 1 1 1 1 1)
COR_CONTIG_LAREPEAT_LEAVE_COV=(1.01 1.01 1.01 1.01 1.01 1.01 1.01)
COR_CONTIG_LAREPEAT_OLEN=0
COR_CONTIG_LAREPEAT_IDENTITYOVLS=1
### TKmerge
COR_CONTIG_TKMERGE_DELETE=1
### TKcombine
COR_CONTIG_TKCOMBINE_DELETE=1
COR_CONTIG_TKCOMBINE_VERBOSE=1
### CTanalyze
COR_CONTIG_CTANALYZE_VERBOSE=1
COR_CONTIG_CTANALYZE_MINCLEN=1000
COR_CONTIG_CTANALYZE_EXPRAWREADCOV=${RAW_COV}
COR_CONTIG_CTANALYZE_MAXSPURLEN=100000
COR_CONTIG_CTANALYZE_MAXTIPLEN=200000
COR_CONTIG_CTANALYZE_REPEATTRACK="repeats_cov${COR_CONTIG_LAREPEAT_COV[3]}_l${COR_CONTIG_LAREPEAT_ENTER_COV[3]}_h${COR_CONTIG_LAREPEAT_LEAVE_COV[3]}_${COR_CONTIG_DATANDER_FOLDER}_dust"
ptype="dalign"
if [[ -n ${FIX_FILT_SCRUB_TYPE} && ${FIX_FILT_SCRUB_TYPE} -eq 2 ]]
then 
    ptype="repcomp"
else 
    ptype="forcealign"
fi
COR_CONTIG_CTANALYZE_TRIMTRACK_PATCHEDREADS="trim1_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}_${ptype}"
COR_CONTIG_CTANALYZE_DIR="analyze01"

# ----------------------------------------------------------------- PACBIO ARROW OPTIONS ----------------------------------------------------------------------------------------------------

### general options
PB_ARROW_RUNID=1                                                          # used for output directory arrow_run${PB_ARROW_RUNID}
PB_ARROW_BAM="/projects/dazzlerAssembly/LAB1608.HYLES_VESPERTILIO/data/pacbio/"    # directory with bam files
PB_ARROW_OUTDIR="${FIX_FILT_OUTDIR}_dalign"
PB_ARROW_INFASTA="${FIX_FILT_OUTDIR}_dalign/correction/contigs"		    # will be ignored if runID is greater then 1
PB_ARROW_MAKEUNIQUEHEADER=0                                                 # to ensure unique header, add sequence index to fasta header 
### pbalign
PB_ARROW_PBALIGN_LOGFILE=1
PB_ARROW_PBALIGN_UNALIGNFILE=1
PB_ARROW_PBALIGN_MINLEN=500
PB_ARROW_PBALIGN_THREADS=24
### arrow
PB_ARROW_ARROW_THREADS=24
PB_ARROW_ARROW_DIPLOID=1
PB_ARROW_ARROW_VERBOSE=1
PB_ARROW_ARROW_REPORTEFFCOV=1
PB_ARROW_ARROW_ANNOTATEGFF=1
PB_ARROW_ARROW_GFFOUT=1
PB_ARROW_ARROW_FQOUT=1
PB_ARROW_ARROW_VCFOUT=1


# ***************************************************************** runtime parameter for slurm settings:  threads, mem, time ***************************************************************

### default parameter for 24-core nodes  #####
THREADS_DEFAULT=1
MEM_DEFAULT=16144
TIME_DEFAULT=24:00:00

##### DBdust  #####
THREADS_DBdust=1
MEM_DBdust=16720
TIME_DBdust=00:30:00
STEPSIZE_DBdust=5

##### datander #####
THREADS_datander=${RAW_REPMASK_DATANDER_THREADS}
MEM_datander=36720
TIME_datander=12:00:00
STEPSIZE_datander=3

#### LAseparate #####
THREADS_LAseparate=1
MEM_LAseparate=16720
TIME_LAseparate=12:30:00
STEPSIZE_LAseparate=100

#### repcomp #####
THREADS_repcomp=${RAW_FIX_REPCOMP_THREADS}
MEM_repcomp=68720
TIME_repcomp=24:00:00
STEPSIZE_repcomp=50

##### forcealign #####
THREADS_forcealign=${RAW_FIX_FORCEALIGN_THREADS}
MEM_forcealign=8720
TIME_forcealign=24:00:00
STEPSIZE_forcealign=10

##### daligner #####
THREADS_daligner=8
MEM_daligner=46720
TIME_daligner=24:00:00

##### LAmerge #####
THREADS_LAmerge=1     
MEM_LAmerge=16144
TIME_LAmerge=24:00:00
TASKS_LAmerge=8

##### LAq #####
THREADS_LAq=1
MEM_LAq=16144
TIME_LAq=12:00:00
TASKS_LAq=12
STEPSIZE_LAq=5

##### TKmerge #####
THREADS_TKmerge=1
MEM_TKmerge=16144
TIME_TKmerge=07:00:00

##### LAstitch #####
THREADS_LAstitch=1
MEM_LAstitch=36144
TIME_LAstitch=12:00:00
TASKS_LAstitch=12

##### LAgap #####
THREADS_LAgap=1
MEM_LAgap=36144
TIME_LAgap=12:00:00
TASKS_LAgap=12

##### TKhomogenize #####
THREADS_TKhomogenize=1
MEM_TKhomogenize=16144
TIME_TKhomogenize=12:00:00
TASKS_TKhomogenize=4
STEPSIZE_TKhomogenize=2

##### OGbuild #####
THREADS_OGbuild=1
MEM_OGbuild=96000
TIME_OGbuild=12:00:00
TASKS_OGbuild=4

##### OGtour #####
THREADS_OGtour=1
MEM_OGtour=64000
TIME_OGtour=12:00:00
TASKS_OGtour=12
STEPSIZE_OGtour=5

##### tour2fasta #####
THREADS_tour2fasta=1
MEM_tour2fasta=64144
TIME_tour2fasta=24:00:00
TASKS_tour2fasta=7
STEPSIZE_tour2fasta=5

##### OGlayout #####
THREADS_OGlayout=1
MEM_OGlayout=16144
TIME_OGlayout=12:00:00
TASKS_OGlayout=24
STEPSIZE_OGlayout=10

##### LAcorrect #####
THREADS_LAcorrect=${FIX_CORR_LACORRECT_THREAD}
MEM_LAcorrect=36000
TIME_LAcorrect=24:00:00

########### pbalign 
THREADS_pbalign=${PB_ARROW_PBALIGN_THREADS}
MEM_pbalign=$((${PB_ARROW_PBALIGN_THREADS}*4096+4096))     ### otherwise bamtools sort will fail !!!!
TIME_pbalign=24:00:00

########### arrow
THREADS_arrow=${PB_ARROW_ARROW_THREADS}
MEM_arrow=$((${PB_ARROW_ARROW_THREADS}*4096+4096))     
TIME_arrow=24:00:00
