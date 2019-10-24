####### MARVEL pipeline variable 

### MARVEL PATH
MARVEL_SOURCE_PATH="/projects/dazzler/pippel/prog/MARVEL/DAMARVEL"
MARVEL_PATH="/projects/dazzler/pippel/prog/MARVEL/DAMARVEL-build/DAmar_TRACE_XOVR_75"
#MARVEL_PATH="/projects/dazzler/pippel/prog/MARVEL/DAMARVEL-build/DAmar_TRACE_XOVR_125"

### REPCOMP PATH
REPCOMP_SOURCE_PATH="/projects/dazzler/pippel/prog/repcomp"
REPCOMP_PATH="LIBMAUS2_DAZZLER_ALIGN_ALIGNMENTFILECONSTANTS_TRACE_XOVR=75 /projects/dazzler/pippel/prog/repcomp-build/repcomp_TRACE_XOVR_75"
#REPCOMP_PATH="LIBMAUS2_DAZZLER_ALIGN_ALIGNMENTFILECONSTANTS_TRACE_XOVR=125 /projects/dazzler/pippel/prog/repcomp-build/repcomp_TRACE_XOVR_125"

### DACCORD PATH - used progs: fastaidrename, forcealign
DACCORD_SOURCE_PATH="/projects/dazzler/pippel/prog/daccord/"
DACCORD_PATH="LIBMAUS2_DAZZLER_ALIGN_ALIGNMENTFILECONSTANTS_TRACE_XOVR=75 /projects/dazzler/pippel/prog/daccord-mpi-build/daccord_0_0_635"
#DACCORD_PATH="LIBMAUS2_DAZZLER_ALIGN_ALIGNMENTFILECONSTANTS_TRACE_XOVR=125 /projects/dazzler/pippel/prog/daccord-mpi-build/daccord_0_0_635"

### LASTOOL PATH - used progs: lassort2
LASTOOLS_SOURCE_PATH="/projects/dazzler/pippel/prog/lastools"
LASTOOLS_PATH="LIBMAUS2_DAZZLER_ALIGN_ALIGNMENTFILECONSTANTS_TRACE_XOVR=75 /projects/dazzler/pippel/prog/lastools-build"
#LASTOOLS_PATH="LIBMAUS2_DAZZLER_ALIGN_ALIGNMENTFILECONSTANTS_TRACE_XOVR=125 /projects/dazzler/pippel/prog/lastools-build"

### DZZLER PATH
DAZZLER_SOURCE_PATH="/projects/dazzler/pippel/prog/dazzlerGIT"
DAZZLER_PATH="/projects/dazzler/pippel/prog/dazzlerGIT/TRACE_XOVR_75"
#DAZZLER_PATH="/projects/dazzler/pippel/prog/dazzlerGIT/TRACE_XOVR_125"

### slurm scripts path
SUBMIT_SCRIPTS_PATH="${MARVEL_PATH}/scripts"

############################## tools for pacbio arrow correction 
CONDA_BASE_ENV="source /projects/dazzler/pippel/prog/miniconda3/bin/activate base"
############################## tools HiC HiGlass pipleine, bwa, samtools, pairstools, cooler, ..;
CONDA_HIC_ENV="source /projects/dazzler/pippel/prog/miniconda3/bin/activate hic"
############################## activate purgehaplotigs environment if requires
CONDA_PURGEHAPLOTIGS_ENV="source /projects/dazzler/pippel/prog/miniconda3/bin/activate purge_haplotigs_env"
############################## tools for whatshap phasing
CONDA_WHATSHAP_ENV="source /projects/dazzler/pippel/prog/miniconda3/bin/activate whatshap"

### ENVIRONMENT VARIABLES 
export PATH=${MARVEL_PATH}/bin:${MARVEL_PATH}/scripts:$PATH
export PYTHONPATH=${MARVEL_PATH}/lib.python:$PYTHONPATH
export SCAFF10X_PATH="/projects/dazzler/pippel/prog/scaffolding/Scaff10X_git/src"
export BIONANO_PATH="/projects/dazzler/pippel/prog/bionano/Solve3.3_10252018"
export SALSA_PATH="/projects/dazzler/pippel/prog/scaffolding/SALSA"
export QUAST_PATH="/projects/dazzler/pippel/prog/quast/"
export JUICER_PATH="/projects/dazzler/pippel/prog/scaffolding/juicer"
export JUICER_TOOLS_PATH="/projects/dazzler/pippel/prog/scaffolding/juicer_tools.1.9.8_jcuda.0.8.jar"
export THREEDDNA_PATH="/projects/dazzler/pippel/prog/scaffolding/3d-dna/"
export LONGRANGER_PATH="/projects/dazzler/pippel/prog/longranger-2.2.2"
export SUPERNOVA_PATH="/projects/dazzler/pippel/prog/supernova-2.1.1"
export ARKS_PATH="/projects/dazzler/pippel/prog/scaffolding/arks-build/bin"
export TIGMINT_PATH="/projects/dazzler/pippel/prog/scaffolding/tigmint/bin"
export LINKS_PATH="/projects/dazzler/pippel/prog/scaffolding/links_v1.8.6/"
export JELLYFISH_PATH="/projects/dazzler/pippel/prog/Jellyfish/jellyfish-2.2.10/bin/"
export GENOMESCOPE_PATH="/projects/dazzler/pippel/prog/genomescope/"
export GATK_PATH="/projects/dazzler/pippel/prog/gatk-4.0.3.0/gatk-package-4.0.3.0-local.jar"
export BCFTOOLS_PATH="/projects/dazzler/pippel/prog/bcftools"
export SEQKIT_PATH="/projects/dazzler/pippel/prog/bin/seqkit"
export FASTP_PATH="/projects/dazzler/pippel/prog/fastp/"

BGZIP_THREADS=6
MARVEL_STATS=1
SLURM_STATS=1

## general information
PROJECT_ID=iHylVes1
GSIZE=600M
SLURM_PARTITION=batch			# default slurm partition - todo define individual partion for tasks
SLURM_NUMACTL=1 

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
TENX_PATH=/projects/dazzlerAssembly/LAB1608.HYLES_VESPERTILIO/data/10x
PATCHING_DIR="patching"
ASSMEBLY_DIR="assembly"
COVERAGE_DIR="coverage"
MITO_DIR="mitochondrion"
QC_DATA_DIR="processedData"

RAW_DALIGN_OUTDIR="dalign"
RAW_REPCOMP_OUTDIR="repcomp"
RAW_DACCORD_OUTDIR="daccord"
RAW_DACCORD_INDIR="dalign"

RAW_REPAMSK_OUTDIR=repmask

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> phase -2 - data QC and statistics and format conversion <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#type-0 [10x - prepare] 						[1-3]: 01_longrangerBasic, 02_longrangerToScaff10Xinput, 03_bxcheck
#type-1 [10x - de novo] 						[1-1]: 01_supernova
#type-2 [10x|HiC - kmer-Gsize estimate] 		[1-2]: 01_genomescope
#type-3 [allData - MASH CONTAMINATION SCREEN] 	[1-5]: 01_mashPrepare, 02_mashSketch, 03_mashCombine, 04_mashPlot, 05_mashScreen
#type-4 [10x - QV]   							[1-4]:  01_QVprepareInput, 02_QVlongrangerAlign, 03_QVcoverage, 04QVqv
RAW_QC_TYPE=0
RAW_QC_SUBMIT_SCRIPTS_FROM=1
RAW_QC_SUBMIT_SCRIPTS_TO=2

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> phase -1 - mitochondrium assembly <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#type-0 steps [1-21]: 1-mitoPrepareInput, 2-mitodaligner, 3-mitoLAmerge, 4-mitoLAfilterMito, 5-mitoPrepareMitoHitDB, 6-mitoHitDBdaligner 7-mitoHitDBLAq 8-mitoHitDBLAfix 09_mitoPrepareMitoHitFixDB 10_mitoHitFixDBdaligner 
#                     11_mitoHitFixDBforcealign 12_mitoHitFixDBLAmerge 13_mitoHitFixDBLAq 14_mitoHitFixDBLAgap 15_mitoHitFixDBLAq 16_mitoHitFixDBLAfilter 17_mitoHitFixDBLAcorrect 18_mitoPrepareMitoHitCorDB,19_mitoHitCorDBdaligner, 20_mitoHitCorDBLAq
#                     21_mitoHitCorDBLAfilter
RAW_MITO_TYPE=0

RAW_MITO_SUBMIT_SCRIPTS_FROM=1
RAW_MITO_SUBMIT_SCRIPTS_TO=21

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> phase 0 - DAScover <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# type-0 steps [1-16]: 1-DBdust, 2-datander, 3-TANmask, 4-CheckTan, 5-rmTan, 6-daligner, 7-CheckDaligner, 8-REPmask, 9-rmDaligner, 10-daligner, 11-CheckDaligner, 12-LAmerge, 13-CheckMerge, 14-rmDaligner, 15-DAScover, 16-REPcover]
RAW_DASCOVER_TYPE=0

RAW_DASCOVER_SUBMIT_SCRIPTS_FROM=1
RAW_DASCOVER_SUBMIT_SCRIPTS_TO=16


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> marvel phase 1 - repeat masking <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# type_0 - stepsp[1-14}: 01_createSubdir, 02_DBdust, 03_Catrack, 04_datander, 05_TANmask, 06_Catrack, 07_daligner, 08_LAmerge, 09_LArepeat, 10_TKmerge, 11-daligner, 12-LAmerge, 13-LArepeat, 14-TKmerge
RAW_REPMASK_TYPE=0

RAW_REPMASK_SUBMIT_SCRIPTS_FROM=1
RAW_REPMASK_SUBMIT_SCRIPTS_TO=9

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> marvel phase 2 - read patching <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 
#type-0 - steps[1-10]: 01-createSubdir, 02-daligner, 03-LAmerge, 04-LArepeat, 05-TKmerge, 06-TKcombine, 07-LAfilter, 08-LAq, 09-TKmerge, 10-LAfix
#type-1 - steps[1-10]: 01-createSubdir, 02-LAseparate, 03-repcomp, 04-LAmerge, 05-LArepeat, 06-TKmerge, 07-TKcombine, 08-LAq, 09-TKmerge, 10-LAfix
#type-2 - steps[1-18]: 01-createSubdir, 02-lassort2, 03-computeIntrinsicQV, 04_Catrack, 05_lasdetectsimplerepeats, 06_mergeAndSortRepeats, 07_lasfilteralignments, 08_mergesym2, 09_filtersym, 10_lasfilteralignmentsborderrepeats, 11_mergesym2, 12_filtersym, 13_filterchainsraw, 14_LAfilterChains, 15_LAfilter, 16_split, 17_LAmerge, 18_LAfix
#type-3 - steps[1-1]:  01_patchStats
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
#type-0 steps [1-3]: 01-createSubdir, 02-LAfilter, 03-LAmerge          ## old pipeline
#type-1 steps [1-15]: 01-createSubdir, 02-lassort2, 03-computeIntrinsicQV, 04_Catrack, 05_lasdetectsimplerepeats, 06_mergeAndSortRepeats, 07_lasfilteralignments, 08_mergesym2, 09_filtersym, 10_lasfilteralignmentsborderrepeats, 11_mergesym2, 12_filtersym, 13_filterchainsraw, 14_LAfilter, 15_LAmerge
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
#type-0 steps: 01_createCorrectedContigDB, 02_DBdust, 03_Catrack, 04_datander, 05_TANmask, 06_Catrack, 07_daligner, 08_LAmerge, 09_LArepeat, 10_TKmerge, 11_TKcombine, 12_LAfilterCTchains, 13_LAmerge, 14_CTanalyze, 15_CTstatistics
COR_CONTIG_SUBMIT_SCRIPTS_FROM=0
COR_CONTIG_SUBMIT_SCRIPTS_TO=15

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> marvel phase 9 - pacbio arrow correction  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

PB_ARROW_TYPE=0
#type-0 steps: 1-prepInFasta, 2-pbalign, 3-bamsplit, 4-bamseparate, 5-bamMerge, 6-arrow, 7-statistics
PB_ARROW_SUBMIT_SCRIPTS_FROM=1
PB_ARROW_SUBMIT_SCRIPTS_TO=7

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> marvel phase 10 - contig purge haplotigs  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

CT_PURGEHAPLOTIGS_TYPE=0
#type-0 steps: 1-prepInFasta, 2-createMinimap2RefIndex, 3-minimap2, 4-bamMerge, 5-readCovHist, 6-contigCovHist, 7-purgeHaplotigs, 8-statistics
CT_PURGEHAPLOTIGS_SUBMIT_SCRIPTS_FROM=1
CT_PURGEHAPLOTIGS_SUBMIT_SCRIPTS_TO=8

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> marvel phase 11 - Freebayes polishing on contigs  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

CT_FREEBAYES_TYPE=1
# Type: 0 [bwa mapping] - 01_FBprepareInput, 02_FBfastp, 03_FBbwa, 04_FBmarkDuplicates, 05_FBfreebayes, 06_FBconsensus, 07_FBstatistics 
# Type: 1 [longranger mapping] - 01_FBprepareInput, 02_FBlongrangerAlign, 03_FBfreebayes, 04_FBconsensus, 05_FBstatistics
CT_FREEBAYES_SUBMIT_SCRIPTS_FROM=1
CT_FREEBAYES_SUBMIT_SCRIPTS_TO=5

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> marvel phase 12 - phasing contigs  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

CT_PHASE_TYPE=1
## type-0 [Whatshap]   - pacbio, 10x: 		01_WhatshapPrepareInput, 02_WhatshapMinimap2PacBio, 03_WhatshapPacBioBamSplitByRef, 04_WhatshapPacBioBamSplitByRef, 05_WhatshapPacBioBamMerge
## type-1 [Longranger] - 10x: 				01_LongrangerPrepareInput, 02_LongrangerLongrangerWgs, 03_LongrangerBcftoolsConsensus, 04_LongrangerStatistics
## type-2 [HapCut2]    - pacbio, 10x, HiC: 	todo
CT_PHASE_SUBMIT_SCRIPTS_FROM=1
CT_PHASE_SUBMIT_SCRIPTS_TO=4

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> marvel phase 13 - scaff10x scaffolding  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

SC_10X_TYPE=2
#type 0: scaff10x - break10x pipeline		steps: 01_scaff10xprepare, 02_scaff10xbreak10, 03_scaff10xscaff10x, 04_scaff10xbreak10x, 05_scaff10xscaff10x, 06_scaff10xbreak10x, 07_scaff10xStatistics
#type 1: tigmint - arks - links pipeline	steps: 01_arksPrepare, 02_arksLongranger, 03_arksTigmint, 04_arksArks, 05_arksLINKS
#type 2: scaff10x using longranger bam		steps: 01_scaff10xprepare, 02_scaff10xLongrangerAlign, 03_scaff10xPrepareIntermediate, 04_scaff10xScaff10x, 05_scaff10xStatistics
#type 3: break10x using longranger bam		steps: 01_break10xPrepare, 02_break10xLongrangerAlign, 03_break10xPrepareIntermediate, 04_break10xBreak10x, 05_break10xStatistics
SC_10X_SUBMIT_SCRIPTS_FROM=1
SC_10X_SUBMIT_SCRIPTS_TO=5

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> marvel phase 14 - bionano scaffolding  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

SC_BIONANO_TYPE=0
# Type: 0 steps: 01_BNscaffold, 02_BNstatistics 
SC_BIONANO_SUBMIT_SCRIPTS_FROM=1
SC_BIONANO_SUBMIT_SCRIPTS_TO=2

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> marvel phase 15 - HiC QC and scaffolding  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

SC_HIC_TYPE=0
# Type: 0 Arima Mapping Pipeline (For QC) 				 steps: 01_HICsalsaPrepareInput, 02_HICsalsaBwa, 03_HICsalsaFilter, 04_HICsalsaMerge, 05_HICsalsaMarkduplicates, 06_HICsalsaSalsa, 07_HICsalsaStatistics 
# Type: 1 Phase Genomics Mapping Pipeline (For QC) 		 steps: 01_HICphasePrepareInput, 02_HICphaseBwa, 03_HICphaseFilter, 04_HICphaseMatlock
# Type: 2 Aiden Lab Juicer/3d-dna Scaffolding Pipeline   steps: 01_HIC3dnaPrepareInput, 02_HIC3dnaJuicer, 03_HIC3dnaAssemblyPipeline
# Type: 3 Aiden Lab Juicer/3d-dna visualization Pipeline steps: 01_HIC3dnaPrepareInput, 02_HIC3dnaJuicer, 03_HIC3dnaVisualize
# Type: 4 - higlass visualization                        steps: 01_HIChiglassPrepare, 02_HiChiglassBwa, 03_HiChiglassFilter, 04_HiChiglassMatrix
SC_HIC_SUBMIT_SCRIPTS_FROM=1
SC_HIC_SUBMIT_SCRIPTS_TO=7
	
# ----------------------------------------------------------------- RAW MITOCHONDRION OPTIONS - always on RAW_DB ---------------------------------------------------------------------------

RAW_MITO_REFFASTA=/projects/dazzlerAssembly/LAB1608.HYLES_VESPERTILIO/data/mitochondria_ref/iHylVes_mt.fasta

# daligner
#RAW_MITO_DALIGNER_KMER=14
RAW_MITO_DALIGNER_ERR=0.7
RAW_MITO_DALIGNER_BIAS=0
RAW_MITO_DALIGNER_OLEN=1000
RAW_MITO_DALIGNER_MEM=16
RAW_MITO_DALIGNER_BLOCKCMP=4
RAW_MITO_DALIGNER_FORBLOCK=1
RAW_MITO_DALIGNER_NUMACTL=1
# LAfilterMito
RAW_MITO_LAFILTERMITO_VERBOSE=0
RAW_MITO_LAFILTERMITO_MINRLEN=4000
RAW_MITO_LAFILTERMITO_MAXRLEN=0			## if unspecified will be set based on reference length minus 1000
RAW_MITO_LAFILTERMITO_UTIPS=1500 		## maximum number of unaligned bases of first and/or last alignment of and alignment chain 
RAW_MITO_LAFILTERMITO_MAXOVH=25			## maximum number of bases, that neighboring alignments of an alignment chain can overlap
RAW_MITO_LAFILTERMITO_MAXGAPLEN=1000	## maximum number of bases, that neighboring alignments of an alignment chain can be separated
RAW_MITO_LAFILTERMITO_PERCCOVLEN=75		## minimum base percentage [50,100], that an alignemnt chain must cover of a read   
#LAq 
RAW_MITO_LAQ_QTRIMCUTOFF=25
RAW_MITO_LAQ_MINSEG=2
### LAfilter
RAW_MITO_LAFILTER_VERBOSE=1
RAW_MITO_LAFILTER_PURGE=1
#RAW_MITO_LAFILTER_OLEN=1000
#RAW_MITO_LAFILTER_RLEN=4000
#RAW_MITO_LAFILTER_DIF=20
RAW_MITO_LAFILTER_UBAS=0
#RAW_MITO_LAFILTER_MINTIPCOV=5
#RAW_MITO_LAFILTER_MULTIMAPPER=2
#RAW_MITO_LAFILTER_REMPERCWORSTALN=20
#RAW_MITO_LAFILTER_EXCLUDEREADS=blacklist.txt
#RAW_MITO_LAFILTER_STITCH=100
#RAW_MITO_LAFILTER_STITCH_AGG=0
RAW_MITO_LAFILTER_TRIM=1
###LAfix
RAW_MITO_LAFIX_GAP=-1
RAW_MITO_LAFIX_MLEN=2000
RAW_MITO_LAFIX_LOW_COVERAGE=0
#RAW_MITO_LAFIX_CONVERTRACKS="track1 track2 ..."
###forcealign
RAW_MITO_FORCEALIGN_PARTIAL=1
RAW_MITO_FORCEALIGN_THREADS=1
#RAW_MITO_FORCEALIGN_MAXDIST=250
#RAW_MITO_FORCEALIGN_BORDER=25
#RAW_MITO_FORCEALIGN_CORRELATION=0.65


# ----------------------------------------------------------------- RAW MASH CONTAMINATION SCREENING ---------------------------------------------------------------------------

MASH_REF_GENOMES=/projects/dazzler/pippel/prog/data/refseq.genomes.k21s1000.msh
RAW_MASH_FASTP_THREADS=1

# ----------------------------------------------------------------- RAW DASCOVER OPTIONS - always on RAW_DAZZ_DB ---------------------------------------------------------------------------

RAW_DASCOVER_DBDUST_BIAS=0

RAW_DASCOVER_DATANDER_THREADS=8
#RAW_DASCOVER_DATANDER_MINLEN=500

DASCOVER_LACHECK_BLOCKCMP=4
DASCOVER_LACHECK_OPT=1
 
DASCOVER_TANMASK_BLOCKCMP=4
DASCOVER_TANMASK_VERBOSE=1

#RAW_DASCOVER_DALIGNER_KMER=14
#RAW_DASCOVER_DALIGNER_ERR=0.7
RAW_DASCOVER_DALIGNER_BIAS=0
#RAW_DASCOVER_DALIGNER_OLEN=1000
RAW_DASCOVER_DALIGNER_MEM=32
DASCOVER_DALIGNER_BLOCKCMP=4
RAW_DASCOVER_DALIGNER_FORBLOCK=1
RAW_DASCOVER_DALIGNER_NUMACTL=1

DASCOVER_REPMASK_BLOCKCMP=4
DASCOVER_REPMASK_VERBOSE=1
DASCOVER_REPMASK_COVERAGE=10

RAW_DASCOVER_LAMERGE_VERBOSE=1

RAW_DASCOVER_DASCOVER_VERBOSE=1
#RAW_DASCOVER_DASCOVER_REPEAT="dust tan rep"
#RAW_DASCOVER_DASCOVER_HGAPLEN=6000

# ----------------------------------------------------------------- RAW REPEAT MASKING OPTIONS - always on RAW_DB ---------------------------------------------------------------------------

# number of daligner block comparison for repeat masking if there are more then one number then mulyiple round of repeat masking are applied
RAW_REPMASK_BLOCKCMP=(1)
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
RAW_FIX_LAREPEAT_MAX_COV=200
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
RAW_FIX_LAFIX_FIXCHIMERS=1
RAW_FIX_LAFIX_DISCARDCHIMERS=1
RAW_FIX_LAFIX_MINCHIMERBORDERCOV=7
RAW_FIX_LAFIX_MAXCHIMERLEN=8000

# ----------------------------------------------------------------- FIX REPEAT MASKING OPTIONS - always on RAW_DB ---------------------------------------------------------------------------

FIX_REPMASK_USELAFIX_PATH=patchedReads_dalign
# number of daligner block comparison for repeat masking if there are more then one number then mulyiple round of repeat masking are applied
FIX_REPMASK_BLOCKCMP=(1)
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
FIX_SCRUB_REPCOMP_MASK="dust"
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
FIX_SCRUB_LAMERGE_SORT=1
FIX_SCRUB_LAMERGE_VERBOSE=0
# LArepeat
FIX_SCRUB_LAREPEAT_LEAVE_COV=(1.5 1.7 1.5 1.7)
FIX_SCRUB_LAREPEAT_ENTER_COV=(2.0 2.0 2.0 2.0)
FIX_SCRUB_LAREPEAT_COV=(-1 -1 ${FIX_COV} ${FIX_COV})
FIX_SCRUB_LAREPEAT_MAX_COV=200
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
FIX_FILT_LAFILTER_MAXREPEATMERGELEN=3000
FIX_FILT_LAFILTER_MAXREPEATMERGEWINDOW=800
#FIX_FILT_LAFILTER_STITCH=${SCRUB_LAGAP_STITCH}
FIX_FILT_LAFILTER_REPEAT_IDX=0
FIX_FILT_LAFILTER_TRIM=1
FIX_FILT_LAFILTER_UBAS=0
FIX_FILT_LAFILTER_PRELOAD=1
FIX_FILT_LAFILTER_DUST="tan_dust"
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
FIX_TOUR_OGTOUR_DEBUG=1
# tour2fasta
FIX_TOUR_2FASTA_TRIM=1
FIX_TOUR_2FASTA_SPLIT=0 ### split at every junction
# OGlayout
FIX_TOUR_OGLAYOUT_VERBOSE=1
FIX_TOUR_OGLAYOUT_DIST=200
FIX_TOUR_OGLAYOUT_RMREVERSEEDGE=1
FIX_TOUR_OGLAYOUT_OUTPUTFORMAT=graphml		### required for Gephi >=0.9, otherwise color attributes are not recognized

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
COR_CONTIG_DALIGNER_KMER=14
COR_CONTIG_DALIGNER_ERR=0.7
COR_CONTIG_DALIGNER_BIAS=0
COR_CONTIG_DALIGNER_RUNID=1
COR_CONTIG_DALIGNER_OLEN=0
COR_CONTIG_DALIGNER_MEM=64
COR_CONTIG_DALIGNER_MASK="dust ${COR_CONTIG_DATANDER_FOLDER}"
COR_CONTIG_DALIGNER_TRACESPACE=0
### LAmerge
COR_CONTIG_LAMERGE_NFILES=255
### LArepeat 
COR_CONTIG_LAREPEAT_COV=(2 3 4 5 50 100 1000)
COR_CONTIG_LAREPEAT_ENTER_COV=(1.01 1.01 1.01 1.01 1.01 1.01 1.01)
COR_CONTIG_LAREPEAT_LEAVE_COV=(0.99 0.99 0.99 0.99 0.99 0.99 0.99)
COR_CONTIG_LAREPEAT_OLEN=0
COR_CONTIG_LAREPEAT_IDENTITYOVLS=1
### LAfilterChains
COR_CONTIG_LAFILTERCHAINS_NREP=500
COR_CONTIG_LAFILTERCHAINS_PURGE=1
COR_CONTIG_LAFILTERCHAINS_KEEP=0              ## keep only best chain
COR_CONTIG_LAFILTERCHAINS_FUZZ=15000
COR_CONTIG_LAFILTERCHAINS_PERCENTCONTAINED=50
COR_CONTIG_LAFILTERCTCHAIN_REPEATTRACK="repeats_c${COR_CONTIG_LAREPEAT_COV[3]}_l${COR_CONTIG_LAREPEAT_LEAVE_COV[3]}h${COR_CONTIG_LAREPEAT_ENTER_COV[3]}_${COR_CONTIG_DATANDER_FOLDER}_dust"
### TKmerge
COR_CONTIG_TKMERGE_DELETE=1
### TKcombine
COR_CONTIG_TKCOMBINE_DELETE=1
COR_CONTIG_TKCOMBINE_VERBOSE=1
### CTanalyze
COR_CONTIG_CTANALYZE_VERBOSE=1
COR_CONTIG_CTANALYZE_EXPRAWREADCOV=${RAW_COV}
COR_CONTIG_CTANALYZE_MAXSPURLEN=200000    
COR_CONTIG_CTANALYZE_MAXTIPLEN=100000        
COR_CONTIG_CTANALYZE_DIR="analyze_01"    	 
COR_CONTIG_CTANALYZE_FUZZYSVLEN=${COR_CONTIG_LAFILTERCHAINS_FUZZ}
COR_CONTIG_CTANALYZE_MINPRIMLEN=150000
COR_CONTIG_CTANALYZE_MINPRIMCREADS=5
COR_CONTIG_CTANALYZE_MAXREPEATPERC=75
ptype="dalign"
if [[ -n ${FIX_FILT_SCRUB_TYPE} && ${FIX_FILT_SCRUB_TYPE} -eq 2 ]]
then 
    ptype="repcomp"
elif [[ -n ${FIX_FILT_SCRUB_TYPE} && ${FIX_FILT_SCRUB_TYPE} -eq 3 ]]
then 
    ptype="forcealign"
fi
COR_CONTIG_CTANALYZE_TRIMTRACK_PATCHEDREADS="trim1_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}_${ptype}"
   
COR_CONTIG_CTANALYZE_CONTIGREPEATTRACK=repeats_c${COR_CONTIG_LAREPEAT_COV[3]}_l${COR_CONTIG_LAREPEAT_LEAVE_COV[3]}h${COR_CONTIG_LAREPEAT_ENTER_COV[3]}_${COR_CONTIG_DATANDER_FOLDER}_dust    
if [[ ${FIX_SCRUB_LAREPEAT_COV[$x]} -ne -1 ]]
then  
      COR_CONTIG_CTANALYZE_READREPEATTRACK=frepeats_c${FIX_SCRUB_LAREPEAT_COV[$x]}_l${FIX_SCRUB_LAREPEAT_LEAVE_COV[$x]}h${FIX_SCRUB_LAREPEAT_ENTER_COV[$x]}_${ptype}_${FIX_REPMASK_LAREPEAT_REPEATTRACK}_${FIX_REPMASK_TANMASK_TRACK}_dust
else  
      COR_CONTIG_CTANALYZE_READREPEATTRACK=frepeats_calCov_l${FIX_SCRUB_LAREPEAT_LEAVE_COV[$x]}h${FIX_SCRUB_LAREPEAT_ENTER_COV[$x]}_${ptype}_${FIX_REPMASK_LAREPEAT_REPEATTRACK}_${FIX_REPMASK_TANMASK_TRACK}_dust
fi  
COR_CONTIG_CTANALYZE_FULLSTATS=0   

# ----------------------------------------------------------------- PACBIO ARROW OPTIONS ----------------------------------------------------------------------------------------------------

### general options
PB_ARROW_RUNID=1                                                          # used for output directory arrow_run${PB_ARROW_RUNID}
PB_ARROW_BAM="${DB_PATH}"   										 # directory with bam files
PB_ARROW_OUTDIR="${FIX_FILT_OUTDIR}"
PB_ARROW_REFFASTA="stats/contigs/m1/haploSplit/test.fa"	
PB_ARROW_P_HEADER="stats/contigs/m1/haploSplit/test.p.header"
PB_ARROW_A_HEADER="stats/contigs/m1/haploSplit/test.a.header"	
PB_ARROW_MAKEUNIQUEHEADER=0                                         # to ensure unique header, add sequence index to fasta header 
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

# ----------------------------------------------------------------- CONTIG FREEBAYES OPTIONS ----------------------------------------------------------------------------------------------------

### general options
CT_FREEBAYES_RUNID=1												# used for output directory purgeHaplotigs_run${PB_ARROW_RUNID}
CT_FREEBAYES_READS="${TENX_PATH}"   								# directory with pacbio fasta files
CT_FREEBAYES_READSTYPE="10x"   								# if its 10x data we need to trim off the first 23 bases
CT_FREEBAYES_OUTDIR="${FIX_FILT_OUTDIR}"
CT_FREEBAYES_REFFASTA="stats/contigs/m1/arrow_2/mMyoMyo_m1_A.p.fasta"	# will be ignored if runID is greater then 1
### fastp
CT_FREEBAYES_FASTP_THREADS=4
### picard tools
CT_FREEBAYES_PICARD_XMX=24						# java memory options in Gb
CT_FREEBAYES_PICARD_XMS=24						# java memory options in Gb
### bwa
CT_FREEBAYES_BWA_THREADS=40
CT_FREEBAYES_BWA_VERBOSITY=3						# 1=error, 2=warning, 3=message, 4+=debugging [3]
### samtools sort
CT_FREEBAYES_SAMTOOLS_THREADS=10
CT_FREEBAYES_SAMTOOLS_MEM=4							# Set maximum memory in Gigabases per thread

# ----------------------------------------------------------------- CONTIG PURGEHAPLPOTIGS OPTIONS ----------------------------------------------------------------------------------------------------

### general options
CT_PURGEHAPLOTIGS_RUNID=1													# used for output directory purgeHaplotigs_run${PB_ARROW_RUNID}
CT_PURGEHAPLOTIGS_PACBIOFASTA="${DB_PATH}"   								# directory with pacbio fasta files
CT_PURGEHAPLOTIGS_OUTDIR="${FIX_FILT_OUTDIR}"
CT_PURGEHAPLOTIGS_INFASTA="stats/contigs/m1/arrow_2/mMyoMyo_m1_A.p.fasta"	# will be ignored if runID is greater then 1
### minimap2
CT_PURGEHAPLOTIGS_MINIMAP2IDXTHREADS=8										# number of threads to create reference index
CT_PURGEHAPLOTIGS_MINIMAP2ALNTHREADS=24										# number of threads to align reads
CT_PURGEHAPLOTIGS_SAMTOOLSTHREADS=8
CT_PURGEHAPLOTIGS_SAMTOOLSMEM=1
### purgeHaplotigs
CT_PURGEHAPLOTIGS_THREADS=24

# ----------------------------------------------------------------- CONTIG WHATSHAP PAHSING OPTIONS ----------------------------------------------------------------------------------------------------

### general whatshap options
CT_PHASE_RUNID=1
CT_PHASE_OUTDIR="${FIX_FILT_OUTDIR}"
CT_PHASE_REFFASTA="stats/contigs/m1/freebayes/mMyoMyo_m1_f.p.fasta"
### use either reads - then full pipeline is started ...
CT_PHASE_READS_10X=${TENX_PATH}
CT_PHASE_READS_PACBIO=${DB_PATH}
CT_PHASE_READS_HIC=${HIC_PATH}
### bwa
CT_PHASE_BWA_THREADS=40
CT_PHASE_BWA_VERBOSITY=3						# 1=error, 2=warning, 3=message, 4+=debugging [3]
### picard tools
CT_PHASE_PICARD_XMX=24						# java memory options in Gb
CT_PHASE_PICARD_XMS=24						# java memory options in Gb
### samtools sort
CT_PHASE_SAMTOOLS_THREADS=10
CT_PHASE_SAMTOOLS_MEM=4							# Set maximum memory in Gigabases per thread
### qv min mapping quality
CT_PHASE_HIC_MINMAPQV=10
CT_PHASE_10X_MINMAPQV=10
CT_PHASE_PACBIO_MINMAPQV=10
 ### minimap2
CT_PHASE_MINIMAP2IDXTHREADS=8										# number of threads to create reference index
CT_PHASE_MINIMAP2ALNTHREADS=24										# number of threads to align reads
CT_PHASE_SAMTOOLSTHREADS=8
CT_PHASE_SAMTOOLSMEM=1

# ----------------------------------------------------------------- SCAFFOLDING - 10X OPTIONS ----------------------------------------------------------------------------------------------------

### general 10x options
SC_10X_RUNID=1
SC_10X_OUTDIR="${FIX_FILT_OUTDIR}"
SC_10X_REF="stats/contigs/m1/freebayes_1/mMyoMyo_m1_f.p.fasta"
SC_10X_READS=${TENX_PATH}
SC_10X_BWA_THREADS=48

### scaff10x options
SC_10X_SCAFF10X_THREADS=48
#SC_10X_SCAFF10X_ALIGNER=bwa		### bwa or smalt
#SC_10X_SCAFF10X_SCORE=20
SC_10X_SCAFF10X_MATRIX=2000
#SC_10X_SCAFF10X_MINREADS=12                            ### VGP: round1: 12, round2: 8 (default: 10)
SC_10X_SCAFF10X_MINREADS_STEP1=12
SC_10X_SCAFF10X_MINREADS_STEP2=8
SC_10X_SCAFF10X_LONGREAD=1
SC_10X_SCAFF10X_GAPSIZE=100                             ### should be the same as used in scaff_reads
SC_10X_SCAFF10X_EDGELEN=50000
#SC_10X_SCAFF10X_MINSHAREDBARCODES=10           ### VGP: round1: 10, round2: 10
SC_10X_SCAFF10X_MINSHAREDBARCODES_STEP1=10
SC_10X_SCAFF10X_MINSHAREDBARCODES_STEP2=10
SC_10X_SCAFF10X_BLOCK=50000				### VGP: round1: 50000, round2: 50000
#SC_10X_SCAFF10X_SAM="path to previously created sam file"
#SC_10X_SCAFF10X_BAM="path to previously created bam file"
#SC_10X_SCAFF10X_READSBC1="m1/scaff10x_1/scaff10x_BC_1.fastq" 		## produced in step1 of scaff10x pipeleine
#SC_10X_SCAFF10X_READSBC2="m1/scaff10x_1/scaff10x_BC_2.fastq"		## produced in step1 of scaff10x pipeleine
### break10x options
SC_10X_BREAK10X_THREADS=48			# nodes  (30)  - number of CPUs requested
SC_10X_BREAK10X_READS=5		       	# reads  (5)   - minimum number of reads per barcode
SC_10X_BREAK10X_SCORE=20      		# score  (20)  - minimum average mapping score on an area covered by reads with the same barcode
SC_10X_BREAK10X_COVER=50       		# cover  (50)  - minimum barcode coverage at the breakpoint
SC_10X_BREAK10X_GAP=100       		# gap    (100) - gap size in building scaffold
SC_10X_BREAK10X_RATIO=15
### tigmint-cut options
SC_10X_TIGMINT_CUT_TRIM=0			#Number of base pairs to trim at contig cuts (bp) [0]
SC_10X_TIGMINT_CUT_SPAN=20			#Spanning molecules threshold (no misassembly in window if num. spanning molecules >= n [2]), but arks uses 20 as default value why??? 
SC_10X_TIGMINT_CUT_WINDOW=1000		#Window size used to check for spanning molecules (bp)  [1000]
### tigmint-molecule options
SC_10X_TIGMINT_MOLECULE_MINMOLSIZE=2000		# Minimum molecule size [2000]
SC_10X_TIGMINT_MOLECULE_ALNSCORERATIO=0.65	#Minimum ratio of alignment score (AS) over read length [0.65]
SC_10X_TIGMINT_MOLECULE_MAXMISMATCH=5		#Maximum number of mismatches (NM) [5]
SC_10X_TIGMINT_MOLECULE_MAXDIST=50000		#Maximum distance between reads in the same molecule [50000]
SC_10X_TIGMINT_MOLECULE_MINMAPQ=0			#Minimum mapping quality [0]
SC_10X_TIGMINT_MOLECULE_NUMREADS=4		#Minimum number of reads per molecule (duplicates are filtered out) [4]  	
### arks
SC_10X_ARKS_RUNTYPE=full				# full, align, graph
SC_10X_ARKS_MINREADPAIRS=5			# Minimum number of mapping read pairs/Index required before creating edge in graph. (default: 5)
SC_10X_ARKS_MULTIPLICITY=50-10000	# Range (in the format min-max) of index multiplicity (only reads with indices in this multiplicity range will be included in graph) (default: 50-10000)
SC_10X_ARKS_MINCONTIGLEN=500		# Minimum contig length to consider for scaffolding (default: 500) 
SC_10X_ARKS_MINFRACTION=0.55		# Minimum fraction of read kmers matching a contigId for a read to be associated with the contigId. (default: 0.55)
SC_10X_ARKS_KMER=30					# k-value for the size of a k-mer. (default: 30) (required)	
SC_10X_ARKS_PVALUE=0.05				# Maximum p-value for H/T assignment and link orientation determination. Lower is more stringent (default: 0.05)	
SC_10X_ARKS_ENDLEN=30000			# End length (bp) of sequences to consider (default: 30000)				
SC_10X_ARKS_DISTESTIMATE=1			# DISTANCE ESTIMATION OPTIONS
SC_10X_ARKS_NUMNEIGHBOUR=20			# num neighbouring samples to estimate distance upper bound [20]			
SC_10X_ARKS_INTRACONTIGTSV=			# -s=FILE output TSV of intra-contig distance/barcode data [disabled]
SC_10X_ARKS_INTERCONTIGTSV=			# -S=FILE output TSV of inter-contig distance/barcode data [disabled]
SC_10X_ARKS_CHECKPOINTS=0			# EXTRA OUTPUT OPTIONS	
									# 0    no checkpoint files (default), 
									# 1    outputs of kmerizing the draft (ContigRecord + ContigKmerMap only)
                        			# 2    output of aligning chromium to draft (IndexMap only)
                        			# 3    all checkpoint files (ContigRecord, ContigKmerMap, and IndexMap)
SC_10X_ARKS_MAXNODEDEGREE=0			# Maximum degree of nodes in graph. All nodes with degree greater than this number will be removed from the graph prior to printing final graph. For no node removal, set to 0 (default: 0)	
SC_10X_ARKS_THREADS=24	
### LINKS
SC_10X_LINKS_VERBOSE=1
SC_10X_LINKS_NUMLINKS=5 			#minimum number of links (k-mer pairs) to compute scaffold (default -l 5, optional)
SC_10X_LINKS_MAXLINKRATIO=0.3		#maximum link ratio between two best contig pairs (default -a 0.3, optional)
SC_10X_LINKS_MINCONTIGLEN=15000		# minimum contig length to consider for scaffolding (default -z 500, optional)
			
# ----------------------------------------------------------------- SCAFFOLDING - BIONANO OPTIONS ----------------------------------------------------------------------------------------------------

### general bionano options
SC_BIONANO_RUNID=1
SC_BIONANO_REF="stats/contigs/m1/scaff10x_1/mMyoMyo_m1_f.p.fasta"
SC_BIONANO_REF_EXCLUDELIST="stats/contigs/m1/haploSplit/filter/mMyoMyo_m1_h.p.excludeP65RepeatContigs.clist"
SC_BIONANO_OUTDIR="${FIX_FILT_OUTDIR}"
SC_BIONANO_FULLSTATS=0

SC_BIONANO_CONFLICTLEVEL_GENOMEMAPS=2			## 1 no filter, 2 cut contig at conflict, 3 exclude conflicting contig
SC_BIONANO_CONFLICTLEVEL_SEQUENCE=2			## 1 no filter, 2 cut contig at conflict, 3 exclude conflicting contig
## single enzyme
SC_BIONANO_ENZYME_1="BSSSI"
SC_BIONANO_ASSEMBLY_1="/projects/dazzlerAssembly/mRhiFer_LAB1237/bionano/Solve3.2.2_withoutReference/mRhiFer_BSSSI_1_fullAssembly_SNRfilt_run01_v322/contigs/mRhiFer_refineFinal1/MRHIFER_REFINEFINAL1.cmap"
# optional parameter, required to map molecules to hybrid scaffold, and compute chimeric quality scores 
SC_BIONANO_MOLCULES_1="/projects/dazzlerAssembly/mRhiFer_LAB1237/bionano/Solve3.2.2_withoutReference/mRhiFer1_Saphyr_BSSSI_SNRfilt.bnx"
SC_BIONANO_ASSEMBLYSCRIPT_1="/projects/dazzlerAssembly/mRhiFer_LAB1237/bionano/Solve3.2.2_withoutReference/clusterArguments_customGPU.xml"
SC_BIONANO_ASSEMBLYOPTARGS_1="/projects/dazzler/pippel/prog/bionano/Solve3.2.2_08222018/RefAligner/7782.7865rel/optArguments_nonhaplotype_saphyr.xml"
SC_BIONANO_ASSEMBLY_NOISE_1="/projects/dazzlerAssembly/mRhiFer_LAB1237/bionano/Solve3.2.2_withoutReference/mRhiFer_BSSSI_1_fullAssembly_SNRfilt_run01_v322/contigs/auto_noise/autoNoise1.errbin"
## two enzyme workflow, optional mapping and chimeric quality score caculation are not supported
#SC_BIONANO_ENZYME_2="BSPQI"
#SC_BIONANO_ASSEMBLY_2="/projects/dazzlerAssembly/mRhiFer_LAB1237/bionano/Solve3.2.2_withoutReference/mRhiFer_BSPQI_1_fullAssembly_SNRfilt_run01_v322/contigs/mRhiFer_refineFinal1/MRHIFER_REFINEFINAL1.cmap"
#SC_BIONANO_CUTS_1=file with valid cut coordinates (absolute file path) 
#SC_BIONANO_CUTS_2=file with valid cut coordinates (absolute file path)
# ----------------------------------------------------------------- SCAFFOLDING - HIC QC AND SALSA, 3DNA, JUICER OPTIONS ----------------------------------------------------------------------------------------------------

### general options
SC_HIC_RUNID=1												# used for output directory purgeHaplotigs_run${PB_ARROW_RUNID}
SC_HIC_READS="${HIC_PATH}"   								# directory with pacbio fasta files
SC_HIC_OUTDIR="${FIX_FILT_OUTDIR}"
SC_HIC_REF="stats/contigs/m1/arrow_2/mMyoMyo_m1_A.p.fasta"	# will be ignored if runID is greater then 1
SC_HIC_REF_EXCLUDELIST="stats/contigs/m1/haploSplit/filter/mMyoMyo_m1_h.p.excludeP65RepeatContigs.clist"
SC_HIC_ENZYME_NAME="Sau3AI"
SC_HIC_ENZYME_SEQ="GATC"
SC_HIC_FULLSTATS=0		
### fastp
SC_HIC_FASTP_THREADS=4
### bwa
SC_HIC_BWA_THREADS=40
SC_HIC_BWA_VERBOSITY=3						# 1=error, 2=warning, 3=message, 4+=debugging [3]
### picard tools
SC_HIC_PICARD_XMX=24						# java memory options in Gb
SC_HIC_PICARD_XMS=24						# java memory options in Gb
### samtools sort
SC_HIC_SAMTOOLS_THREADS=10
SC_HIC_SAMTOOLS_MEM=4							# Set maximum memory in Gigabases per thread
### arima qv min mapping quality
SC_HIC_MINMAPQV=10
### juicer and 3d-dna HiC options 
#SC_HIC_JUICER_STAGE=		# Can be: [merge, dedup, final, postproc, early]
SC_HIC_JUICER_SHORTQUEUE=batch
SC_HIC_JUICER_SHORTQUEUETLIMIIT=1200
SC_HIC_JUICER_LONGQUEUE=long
SC_HIC_JUICER_LONGQUEUETLIMIIT=3600
SC_HIC_JUICER_CHUNKSIZE=60000000		#number of lines in split files, must be multiple of 4 (default 90000000, which equals 22.5 million reads
SC_HIC_3DDNA_MODE=haploid 				#Can be: [haploid, diploid]
SC_HIC_3DDNA_MINCONTIGLEN=15000			#Specifies threshold input contig/scaffold size (default is 15000). Contigs/scaffolds smaller than input_size are going to be ignored.
SC_HIC_3DDNA_ROUNDS=2					#Specifies number of iterative rounds for misjoin correction (default is 2).
#SC_HIC_3DDNA_STAGE=					#can be: [polish, split, seal, merge, finalize
#SC_HIC_3DDNA_MAPQV=1					#Mapq threshold for scaffolding and visualization (default is 1).
### 3d-dna visualize pipeline
SC_HIC_3DDNAVISUALIZE_MAPQV=1			#Build map for a specific mapq threshold (default is 1).
#SC_HIC_3DDNAVISUALIZE_MNDPATH=""		#Path to mnd already lifted from input to assembly chromosome: used to skip the remapping step.	
#SC_HIC_3DDNAVISUALIZE_SKIPNORM=1		#Skip normalization.
#SC_HIC_3DDNAVISUALIZE_RESOLUTION= 		#Build for specific resolutions (default is -r 2500000,1000000,500000,250000,100000,50000,25000,10000,5000,1000)
SC_HIC_3DDNAVISUALIZE_CLEANUP=1			#Clean up when done (default: no cleanup.)
SC_HIC_3DDNAVISUALIZE_IGNOREMAPQV=0		#Ignore mapq suffix.
### HiGlass pipeline
SC_HIC_HIGLASS_COOLERRESOLUTION=(1000 5000)	# cooler binning: binsize : e.g.) 5000 (high resolution), 500000 (lower resolution)
SC_HIC_HIGLASS_PAIRTOOLSTHREADS=24
SC_HIC_FULLSTATS=0
### SALSA
SC_HIC_REF_DOBREAKS=1
SC_HIC_REF_NUMITER=3       	
SC_HIC_MINCONTIG=1000 
SC_HIC_KEEPINTERMEDIATES=1       	

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

########### minimap2 alignment
THREADS_minimap2=${CT_PURGEHAPLOTIGS_MINIMAP2ALNTHREADS}
MEM_minimap2=$((${CT_PURGEHAPLOTIGS_MINIMAP2ALNTHREADS}*4096))     
TIME_minimap2=24:00:00

########### minimap2 reference index
THREADS_createMinimap2RefIndex=${CT_PURGEHAPLOTIGS_MINIMAP2IDXTHREADS}
MEM_createMinimap2RefIndex=$((${CT_PURGEHAPLOTIGS_MINIMAP2IDXTHREADS}*4096))     
TIME_createMinimap2RefIndex=24:00:00

########### samtools merge
THREADS_bamMerge=${CT_PURGEHAPLOTIGS_SAMTOOLSTHREADS}
MEM_bamMerge=$((${CT_PURGEHAPLOTIGS_SAMTOOLSTHREADS}*4096))     
TIME_bamMerge=24:00:00

########### purgeHaplotigs readCovHist
THREADS_readCovHist=${CT_PURGEHAPLOTIGS_THREADS}
MEM_readCovHist=$((${CT_PURGEHAPLOTIGS_THREADS}*2048))     
TIME_readCovHist=24:00:00

########### freebayes fastp 
THREADS_FBfastp=${CT_FREEBAYES_FASTP_THREADS}
MEM_FBfastp=$((${CT_FREEBAYES_FASTP_THREADS}*4096))     
TIME_FBfastp=24:00:00

########### freebayes bwa 
THREADS_FBbwa=${CT_FREEBAYES_BWA_THREADS}
MEM_FBbwa=$((${CT_FREEBAYES_BWA_THREADS}*1024+${CT_FREEBAYES_SAMTOOLS_THREADS}*${CT_FREEBAYES_SAMTOOLS_MEM}*1024))     
TIME_FBbwa=24:00:00

########### freebayes picard markduplicates - run indexing parallel  
THREADS_FBmarkDuplicates=${CT_FREEBAYES_SAMTOOLS_THREADS}
MEM_FBmarkDuplicates=$((${CT_FREEBAYES_SAMTOOLS_THREADS}*${CT_FREEBAYES_SAMTOOLS_MEM}*1024+${CT_FREEBAYES_PICARD_XMS}*1024))     
TIME_FBmarkDuplicates=24:00:00

########### HIC salsa bwa 
THREADS_HICsalsaBwa=${SC_HIC_BWA_THREADS}
MEM_HICsalsaBwa=$((${SC_HIC_BWA_THREADS}*1024+${SC_HIC_SAMTOOLS_THREADS}*${SC_HIC_SAMTOOLS_MEM}*1024))
TIME_HICsalsaBwa=24:00:00

########### HIC picard markduplicates - run indexing parallel  
THREADS_HICsalsaMarkduplicates=${SC_HIC_SAMTOOLS_THREADS}
MEM_HICsalsaMarkduplicates=$((${SC_HIC_SAMTOOLS_THREADS}*${SC_HIC_SAMTOOLS_MEM}*1024+${SC_HIC_PICARD_XMS}*1024))
TIME_HICsalsaMarkduplicates=24:00:00

########### WHATSHAP pipeline
THREADS_WHprepareInput=1
MEM_WHprepareInput=46000	     
TIME_WHprepareInput=24:00:00

THREADS_WHPacBioMinimap2=${CT_PHASE_MINIMAP2ALNTHREADS}
MEM_WHPacBioMinimap2=$((${CT_PHASE_MINIMAP2ALNTHREADS}*4096))     
TIME_WHPacBioMinimap2=24:00:00

##### MITO PIPELINE
THREADS_mitodaligner=4
MEM_mitodaligner=$((24*1024))
TIME_mitodaligner=04:00:00

##### SCAFF10X pipeline
THREADS_scaff10Xscaff10x=${SC_10X_SCAFF10X_THREADS}
MEM_scaff10Xscaff10x=$((${SC_10X_SCAFF10X_THREADS}*4096))
TIME_scaff10Xscaff10x=24:00:00

THREADS_scaff10Xbreak10x=${SC_10X_BREAK10X_THREADS}
MEM_scaff10Xbreak10x=$((${SC_10X_BREAK10X_THREADS}*4096))
TIME_scaff10Xbreak10x=24:00:00

##### Bionano pipeline
THREADS_BNscaffold=24
if [[ "${SLURM_PARTITION}" == "gpu" ]]
then 
	THREADS_BNscaffold=40
elif [[ "${SLURM_PARTITION}" == "bigmem" ]]
then 
	THREADS_BNscaffold=48
fi
MEM_BNscaffold=96000
TIME_BNscaffold=24:00:00

#### Juicer/3d-dna pipeline
THREADS_juicer=24
if [[ "${SLURM_PARTITION}" == "gpu" ]]
then 
	THREADS_juicer=40
elif [[ "${SLURM_PARTITION}" == "bigmem" ]]
then 
	THREADS_juicer=48
fi

THREADS_HIC3dnaJuicer=${THREADS_juicer}
MEM_HIC3dnaJuicer=$((${THREADS_juicer}*4096))
TIME_HIC3dnaJuicer=24:00:00

THREADS_HIC3dnaAssemblyPipeline=${THREADS_juicer}
MEM_HIC3dnaAssemblyPipeline=$((${THREADS_juicer}*4096))
TIME_HIC3dnaAssemblyPipeline=24:00:00

THREADS_HIC3dnaVisualizePipeline=${THREADS_juicer}
MEM_HIC3dnaVisualizePipeline=$((${THREADS_juicer}*4096))
TIME_HIC3dnaVisualizePipeline=24:00:00

THREADS_mashPlot=1
MEM_mashPlot=64000
TIME_mashPlot=24:00:00

THREADS_mashScreen=${THREADS_juicer}
MEM_mashScreen=64000
TIME_mashScreen=24:00:00

THREADS_arksTigmint=${SC_10X_BWA_THREADS}
MEM_arksTigmint=$((${SC_10X_BWA_THREADS}*4096))
TIME_arksTigmint=24:00:00

THREADS_longrangerBasic=${THREADS_juicer}
MEM_longrangerBasic=$((${THREADS_juicer}*8192))
TIME_longrangerBasic=24:00:00

#THREADS_FBlongrangerAlign=${THREADS_juicer}
#MEM_FBlongrangerAlign=$((${THREADS_juicer}*8192))
#TIME_FBlongrangerAlign=240:00:00

THREADS_FBfreebayes=1
MEM_FBfreebayes=32768
TIME_FBfreebayes=24:00:00

THREADS_supernova=${THREADS_juicer}
MEM_supernova=$((${THREADS_juicer}*8192))
TIME_supernova=24:00:00

## HiGlass pipeline
THREADS_HiChiglassFilter=${SC_HIC_HIGLASS_PAIRTOOLSTHREADS}
MEM_HiChiglassFilter=64000
TIME_HiChiglassFilter=24:00:00

THREADS_HiChiglassBwa=${SC_HIC_BWA_THREADS}
MEM_HiChiglassBwa=64000
TIME_HiChiglassBwa=24:00:00

THREADS_HiChiglassMatrix=${SC_HIC_HIGLASS_PAIRTOOLSTHREADS}
MEM_HiChiglassMatrix=64000
TIME_HiChiglassMatrixs=24:00:00

########## statistics
if [[ -z ${BGZIP_THREADS} ]]
then
	BGZIP_THREADS=1
fi
THREADS_statistics=${BGZIP_THREADS}
MEM_statistics=$((${BGZIP_THREADS}*4096))     
TIME_statistics=24:00:00
