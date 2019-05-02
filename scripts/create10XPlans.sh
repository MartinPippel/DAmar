#!/bin/bash -e

configFile=$1
currentStep=$2
slurmID=$3

if [[ ! -f ${configFile} ]]
then 
    (>&2 echo "cannot access config file ${configFile}")
    exit 1
fi

source ${configFile}

gsize=${GSIZE}
i=$((${#GSIZE}-1))
if [[ "${GSIZE: -1}" =~ [gG] ]]
then
 gsize=$((${GSIZE:0:$i}*1000*1000*1000))
fi
if [[ "${GSIZE: -1}" =~ [mM] ]]
then
 gsize=$((${GSIZE:0:$i}*1000*1000))
fi
if [[ "${GSIZE: -1}" =~ [kK] ]]
then
 gsize=$((${GSIZE:0:$i}*1000))
fi

if [[ -z ${TIGMINT_PATH} || ! -f ${TIGMINT_PATH}/tigmint-cut ]]
then
	(>&2 echo "Variable TIGMINT_PATH must be set to a proper tigmint installation directory!!")
    exit 1
fi

if [[ -z ${SCAFF10X_PATH} || ! -f ${SCAFF10X_PATH}/scaff_reads ]]
then
	(>&2 echo "Variable SCAFF10X_PATH must be set to a proper scaff10x installation directory!!")
    exit 1
fi

if [[ -z ${QUAST_PATH} || ! -f ${QUAST_PATH}/quast.py ]]
then
	(>&2 echo "Variable QUAST_PATH must be set to a proper quast installation directory!!")
    exit 1
fi

if [[ -z ${LONGRANGER_PATH} || ! -f ${LONGRANGER_PATH}/longranger ]]
then
	(>&2 echo "Variable LONGRANGER_PATH must be set to a proper longranger installation directory!!")
    exit 1
fi

if [[ -z ${ARKS_PATH} || ! -f ${ARKS_PATH}/arks ]]
then
	(>&2 echo "Variable ARKS_PATH must be set to a proper ARKS installation directory!!")
    exit 1
fi

if [[ -z ${LINKS_PATH} || ! -f ${LINKS_PATH}/LINKS ]]
then
	(>&2 echo "Variable LINKS_PATH must be set to a proper LINKS installation directory!!")
    exit 1
fi

function setScaff10xOptions()
{
	SCAFF10X_SCAFF10X_OPT=""
	
	if [[ -n ${SC_10X_SCAFF10X_THREADS} && ${SC_10X_SCAFF10X_THREADS} -ne 0 ]]
	then
		SCAFF10X_SCAFF10X_OPT="${SCAFF10X_SCAFF10X_OPT} -nodes ${SC_10X_SCAFF10X_THREADS}"		
	fi
	
	if [[ -n ${SC_10X_SCAFF10X_ALIGNER} ]]
	then
		SCAFF10X_SCAFF10X_OPT="${SCAFF10X_SCAFF10X_OPT} -align ${SC_10X_SCAFF10X_ALIGNER}"		
	fi
	
	if [[ -n ${SC_10X_SCAFF10X_SCORE} && ${SC_10X_SCAFF10X_SCORE} -ne 0 ]]
	then
		SCAFF10X_SCAFF10X_OPT="${SCAFF10X_SCAFF10X_OPT} -score ${SC_10X_SCAFF10X_SCORE}"		
	fi
	
	if [[ -n ${SC_10X_SCAFF10X_MATRIX} && ${SC_10X_SCAFF10X_MATRIX} -ne 0 ]]
	then
		SCAFF10X_SCAFF10X_OPT="${SCAFF10X_SCAFF10X_OPT} -matrix ${SC_10X_SCAFF10X_MATRIX}"		
	fi
	
	if [[ -n ${SC_10X_SCAFF10X_MINREADS} && ${SC_10X_SCAFF10X_MINREADS} -ne 0 ]]
	then
		SCAFF10X_SCAFF10X_OPT="${SCAFF10X_SCAFF10X_OPT} -reads ${SC_10X_SCAFF10X_MINREADS}"		
	fi
	
	if [[ -n ${SC_10X_SCAFF10X_MINREADS_STEP1} && ${SC_10X_SCAFF10X_MINREADS_STEP1} -ne 0 ]]
	then
		SCAFF10X_SCAFF10X_OPT="${SCAFF10X_SCAFF10X_OPT} -read-s1 ${SC_10X_SCAFF10X_MINREADS_STEP1}"		
	fi
	
	if [[ -n ${SC_10X_SCAFF10X_MINREADS_STEP2} && ${SC_10X_SCAFF10X_MINREADS_STEP2} -ne 0 ]]
	then
		SCAFF10X_SCAFF10X_OPT="${SCAFF10X_SCAFF10X_OPT} -read-s2 ${SC_10X_SCAFF10X_MINREADS_STEP2}"		
	fi
	
	if [[ -n ${SC_10X_SCAFF10X_LONGREAD} && ${SC_10X_SCAFF10X_LONGREAD} -ne 0 ]]
	then
		SCAFF10X_SCAFF10X_OPT="${SCAFF10X_SCAFF10X_OPT} -longread ${SC_10X_SCAFF10X_LONGREAD}"		
	fi
	
	if [[ -n ${SC_10X_SCAFF10X_GAPSIZE} && ${SC_10X_SCAFF10X_GAPSIZE} -ne 0 ]]
	then
		SCAFF10X_SCAFF10X_OPT="${SCAFF10X_SCAFF10X_OPT} -gap ${SC_10X_SCAFF10X_GAPSIZE}"		
	fi
	
	if [[ -n ${SC_10X_SCAFF10X_EDGELEN} && ${SC_10X_SCAFF10X_EDGELEN} -ne 0 ]]
	then
		SCAFF10X_SCAFF10X_OPT="${SCAFF10X_SCAFF10X_OPT} -edge ${SC_10X_SCAFF10X_EDGELEN}"		
	fi
	
	if [[ -n ${SC_10X_SCAFF10X_MINSHAREDBARCODES} && ${SC_10X_SCAFF10X_MINSHAREDBARCODES} -ne 0 ]]
	then
		SCAFF10X_SCAFF10X_OPT="${SCAFF10X_SCAFF10X_OPT} -link ${SC_10X_SCAFF10X_MINSHAREDBARCODES}"		
	fi

	if [[ -n ${SC_10X_SCAFF10X_MINSHAREDBARCODES_STEP1} && ${SC_10X_SCAFF10X_MINSHAREDBARCODES_STEP1} -ne 0 ]]
	then
		SCAFF10X_SCAFF10X_OPT="${SCAFF10X_SCAFF10X_OPT} -link-s1 ${SC_10X_SCAFF10X_MINSHAREDBARCODES_STEP1}"		
	fi

	if [[ -n ${SC_10X_SCAFF10X_MINSHAREDBARCODES_STEP2} && ${SC_10X_SCAFF10X_MINSHAREDBARCODES_STEP2} -ne 0 ]]
	then
		SCAFF10X_SCAFF10X_OPT="${SCAFF10X_SCAFF10X_OPT} -link-s2 ${SC_10X_SCAFF10X_MINSHAREDBARCODES_STEP2}"		
	fi
	
	if [[ -n ${SC_10X_SCAFF10X_BLOCK} && ${SC_10X_SCAFF10X_BLOCK} -ne 0 ]]
	then
		SCAFF10X_SCAFF10X_OPT="${SCAFF10X_SCAFF10X_OPT} -block ${SC_10X_SCAFF10X_BLOCK}"		
	fi
	
    ### check input variable variables, and overwrite default pipeline if required
    if [[ -n ${SC_10X_SCAFF10X_SAM} && -f ${SC_10X_SCAFF10X_SAM} ]]
    then
    	SCAFF10X_SCAFF10X_OPT="${SCAFF10X_SCAFF10X_OPT} -sam ${SC_10X_SCAFF10X_SAM}"
    elif [[ -n ${SC_10X_SCAFF10X_BAM} && -f ${SC_10X_SCAFF10X_BAM} ]]
    then
    	SCAFF10X_SCAFF10X_OPT="${SCAFF10X_SCAFF10X_OPT} -bam ${SC_10X_SCAFF10X_BAM}"	
	fi
}

function setBreak10xOptions()
{
	SCAFF10X_BREAK10X_OPT=""
	
	if [[ -n ${SC_10X_BREAK10X_THREADS} && ${SC_10X_BREAK10X_THREADS} -ne 0 ]]
	then
		SCAFF10X_BREAK10X_OPT="${SCAFF10X_BREAK10X_OPT} -nodes ${SC_10X_BREAK10X_THREADS}"		
	fi
	
	if [[ -n ${SC_10X_BREAK10X_READS} && ${SC_10X_BREAK10X_READS} -ne 0 ]]
	then
		SCAFF10X_BREAK10X_OPT="${SCAFF10X_BREAK10X_OPT} -reads ${SC_10X_BREAK10X_READS}"		
	fi

	if [[ -n ${SC_10X_BREAK10X_SCORE} && ${SC_10X_BREAK10X_SCORE} -ne 0 ]]
	then
		SCAFF10X_BREAK10X_OPT="${SCAFF10X_BREAK10X_OPT} -score ${SC_10X_BREAK10X_SCORE}"		
	fi

	if [[ -n ${SC_10X_BREAK10X_COVER} && ${SC_10X_BREAK10X_COVER} -ne 0 ]]
	then
		SCAFF10X_BREAK10X_OPT="${SCAFF10X_BREAK10X_OPT} -cover ${SC_10X_BREAK10X_COVER}"		
	fi
	
	if [[ -n ${SC_10X_BREAK10X_RATIO} && ${SC_10X_BREAK10X_RATIO} -ne 0 ]]
	then
		SCAFF10X_BREAK10X_OPT="${SCAFF10X_BREAK10X_OPT} -ratio ${SC_10X_BREAK10X_RATIO}"		
	fi	
	
	if [[ -n ${SC_10X_BREAK10X_GAP} && ${SC_10X_BREAK10X_GAP} -ne 0 ]]
	then
		SCAFF10X_BREAK10X_OPT="${SCAFF10X_BREAK10X_OPT} -gap ${SC_10X_BREAK10X_GAP}"		
	fi
}

function setbwaOptions()
{
	# by default use -p smart pairing, i.e. only one input file and -C append FASTA/FASTQ comment to SAM output, i.e. append 10x barcode
	SCAFFOLD_BWA_OPT=" -p -C"
	
	if [[ -z "${SC_10X_BWA_THREADS}" ]]
	then 
		SC_10X_BWA_THREADS=4	
	fi
	
	SCAFFOLD_BWA_OPT="${SCAFFOLD_BWA_OPT} -t ${SC_10X_BWA_THREADS}"
	
	if [[ -n ${SC_10X_BWA_VERBOSITY} ]]
	then 
		SCAFFOLD_BWA_OPT="${SCAFFOLD_BWA_OPT} -v ${SC_10X_BWA_VERBOSITY}"
	fi
	
	if [[ -n ${SC_10X_BWA_MISMATCHPENALTY} && ${SC_10X_BWA_MISMATCHPENALTY} -gt 0 ]]
	then 
		SCAFFOLD_BWA_OPT="${SCAFFOLD_BWA_OPT} -B ${SC_10X_BWA_MISMATCHPENALTY}"
	fi
}

function setSamtoolsOptions()
{
	SCAFFOLD_SAMTOOLS_OPT=" -tBX"
	
	if [[ -z "${SC_10X_SAMTOOLS_THREADS}" ]]
	then 
		SC_10X_SAMTOOLS_THREADS=4	
	fi
	
	SCAFFOLD_SAMTOOLS_OPT="${SCAFFOLD_SAMTOOLS_OPT} -@ ${SC_10X_SAMTOOLS_THREADS}"
	
	if [[ -n ${SC_10X_SAMTOOLS_MEM} ]]
	then 
		SCAFFOLD_SAMTOOLS_OPT="${SCAFFOLD_SAMTOOLS_OPT} -m ${SC_10X_SAMTOOLS_MEM}G"
	else
		SC_10X_SAMTOOLS_MEM=2
		SCAFFOLD_SAMTOOLS_OPT="${SCAFFOLD_SAMTOOLS_OPT} -m ${SC_10X_SAMTOOLS_MEM}G"
	fi		
}

function setTigmintOptions()
{
	TIGMINT_MOLECULE_OPT=""
	TIGMINT_CUT_OPT=""
	
	# ensure that bwa threads are set
	if [[ -z "${SC_10X_BWA_THREADS}" ]]
	then 
		setbwaOptions
	fi
	
	### set tigmint-cut options
	TIGMINT_CUT_OPT="${TIGMINT_CUT_OPT} --processes ${SC_10X_BWA_THREADS}"
	#Number of base pairs to trim at contig cuts (bp) [0]
	if [[ -n ${SC_10X_TIGMINT_CUT_TRIM} && ${SC_10X_TIGMINT_CUT_TRIM} -gt 0 ]]
	then
		TIGMINT_CUT_OPT="${TIGMINT_CUT_OPT} --trim ${SC_10X_TIGMINT_CUT_TRIM}"	
	fi
	#Spanning molecules threshold (no misassembly in window if num. spanning molecules >= n [2]), but arks uses 20 as default value why??? 
	if [[ -n ${SC_10X_TIGMINT_CUT_SPAN} && ${SC_10X_TIGMINT_CUT_SPAN} -gt 0 ]]
	then
		TIGMINT_CUT_OPT="${TIGMINT_CUT_OPT} --spanning ${SC_10X_TIGMINT_CUT_SPAN}"	
	fi
	#Window size used to check for spanning molecules (bp)  [1000]
	if [[ -n ${SC_10X_TIGMINT_CUT_WINDOW} && ${SC_10X_TIGMINT_CUT_WINDOW} -gt 0 ]]
	then
		TIGMINT_CUT_OPT="${TIGMINT_CUT_OPT} --window ${SC_10X_TIGMINT_CUT_WINDOW}"	
	fi
		
	
	### set tigmint-molecule options
	# Minimum molecule size [2000]
	if [[ -n ${SC_10X_TIGMINT_MOLECULE_MINMOLSIZE} && ${SC_10X_TIGMINT_MOLECULE_MINMOLSIZE} -gt 0 ]]
	then
		TIGMINT_MOLECULE_OPT="${TIGMINT_MOLECULE_OPT} --size ${SC_10X_TIGMINT_MOLECULE_MINMOLSIZE}"
	else
		SC_10X_TIGMINT_MOLECULE_MINMOLSIZE=2000 	
		TIGMINT_MOLECULE_OPT="${TIGMINT_MOLECULE_OPT} --size ${SC_10X_TIGMINT_MOLECULE_MINMOLSIZE}"
	fi
	#Minimum ratio of alignment score (AS) over read length [0.65]
	if [[ -n ${SC_10X_TIGMINT_MOLECULE_ALNSCORERATIO} ]]
	then
		TIGMINT_MOLECULE_OPT="${TIGMINT_MOLECULE_OPT} --as-ratio ${SC_10X_TIGMINT_MOLECULE_ALNSCORERATIO}"
	else
		SC_10X_TIGMINT_MOLECULE_ALNSCORERATIO=0.65	
		TIGMINT_MOLECULE_OPT="${TIGMINT_MOLECULE_OPT} --as-ratio ${SC_10X_TIGMINT_MOLECULE_ALNSCORERATIO}"
	fi
	#Maximum number of mismatches (NM) [5]
	if [[ -n ${SC_10X_TIGMINT_MOLECULE_MAXMISMATCH} && ${SC_10X_TIGMINT_MOLECULE_MAXMISMATCH} -gt 0 ]]
	then
		TIGMINT_MOLECULE_OPT="${TIGMINT_MOLECULE_OPT} --nm ${SC_10X_TIGMINT_MOLECULE_MAXMISMATCH}"
	else
		SC_10X_TIGMINT_MOLECULE_MAXMISMATCH=5
		TIGMINT_MOLECULE_OPT="${TIGMINT_MOLECULE_OPT} --nm ${SC_10X_TIGMINT_MOLECULE_MAXMISMATCH}"
	fi    
	#Maximum distance between reads in the same molecule [50000]
	if [[ -n ${SC_10X_TIGMINT_MOLECULE_MAXDIST} && ${SC_10X_TIGMINT_MOLECULE_MAXDIST} -gt 0 ]]
	then
		TIGMINT_MOLECULE_OPT="${TIGMINT_MOLECULE_OPT} --dist ${SC_10X_TIGMINT_MOLECULE_MAXDIST}"	
	fi   
	#Minimum mapping quality [0]
    if [[ -n ${SC_10X_TIGMINT_MOLECULE_MINMAPQ} && ${SC_10X_TIGMINT_MOLECULE_MINMAPQ} -gt 0 ]]
	then
		TIGMINT_MOLECULE_OPT="${TIGMINT_MOLECULE_OPT} --mapq ${SC_10X_TIGMINT_MOLECULE_MINMAPQ}"	
	fi
	#Minimum number of reads per molecule (duplicates are filtered out) [4]  		
	if [[ -n ${SC_10X_TIGMINT_MOLECULE_NUMREADS} && ${SC_10X_TIGMINT_MOLECULE_NUMREADS} -gt 0 ]]
	then
		TIGMINT_MOLECULE_OPT="${TIGMINT_MOLECULE_OPT} --reads ${SC_10X_TIGMINT_MOLECULE_NUMREADS}"	
	fi											
}

function setArksOptions()
{
	
	ARKS_OPT=""
	
	#Minimum number of mapping read pairs/Index required before creating edge in graph. (default: 5)
	if [[ -n ${SC_10X_ARKS_MINREADPAIRS} && ${SC_10X_ARKS_MINREADPAIRS} -gt 0 ]]
	then
		ARKS_OPT="${ARKS_OPT} -c ${SC_10X_ARKS_MINREADPAIRS}"
	else
		SC_10X_ARKS_MINREADPAIRS=5
		ARKS_OPT="${ARKS_OPT} -c ${SC_10X_ARKS_MINREADPAIRS}" 
	fi
	
	#Range (in the format min-max) of index multiplicity (only reads with indices in this multiplicity range will be included in graph) (default: 50-10000)
	if [[ -n ${SC_10X_ARKS_MULTIPLICITY} ]]
	then
		ARKS_OPT="${ARKS_OPT} -m ${SC_10X_ARKS_MULTIPLICITY}"
	else
		SC_10X_ARKS_MULTIPLICITY="50-10000"
		ARKS_OPT="${ARKS_OPT} -m ${SC_10X_ARKS_MULTIPLICITY}"
	fi

	#Minimum contig length to consider for scaffolding (default: 500)
	if [[ -n ${SC_10X_ARKS_MINCONTIGLEN} && ${SC_10X_ARKS_MINCONTIGLEN} -gt 0 ]]
	then
		ARKS_OPT="${ARKS_OPT} -z ${SC_10X_ARKS_MINCONTIGLEN}"
	else
		SC_10X_ARKS_MINCONTIGLEN=500
		ARKS_OPT="${ARKS_OPT} -z ${SC_10X_ARKS_MINCONTIGLEN}"
	fi
	
	#Minimum fraction of read kmers matching a contigId for a read to be associated with the contigId. (default: 0.55)
	if [[ -n ${SC_10X_ARKS_MINFRACTION} ]]
	then
		ARKS_OPT="${ARKS_OPT} -j ${SC_10X_ARKS_MINFRACTION}"
	fi
	
	# k-value for the size of a k-mer. (default: 30) (required)
	if [[ -n ${SC_10X_ARKS_KMER} && ${SC_10X_ARKS_KMER} -gt 0 ]]
	then
		ARKS_OPT="${ARKS_OPT} -k ${SC_10X_ARKS_KMER}"
	else
		SC_10X_ARKS_KMER=30
		ARKS_OPT="${ARKS_OPT} -k ${SC_10X_ARKS_KMER}"		
	fi

	#Maximum p-value for H/T assignment and link orientation determination. Lower is more stringent (default: 0.05)
	if [[ -n ${SC_10X_ARKS_PVALUE} ]]
	then
		ARKS_OPT="${ARKS_OPT} -r ${SC_10X_ARKS_PVALUE}"
	else
		SC_10X_ARKS_PVALUE=0.05
		ARKS_OPT="${ARKS_OPT} -r ${SC_10X_ARKS_PVALUE}"
	fi	
	
	# End length (bp) of sequences to consider (default: 30000)
	if [[ -n ${SC_10X_ARKS_ENDLEN} && ${SC_10X_ARKS_ENDLEN} -gt 0 ]]
	then
		ARKS_OPT="${ARKS_OPT} -e ${SC_10X_ARKS_ENDLEN}"
	else
		SC_10X_ARKS_ENDLEN=30000
		ARKS_OPT="${ARKS_OPT} -e ${SC_10X_ARKS_ENDLEN}"
	fi
	
	# DISTANCE ESTIMATION OPTIONS
	if [[ -n ${SC_10X_ARKS_DISTESTIMATE} && ${SC_10X_ARKS_DISTESTIMATE} -gt 0 ]]
	then
		ARKS_OPT="${ARKS_OPT} -D"
		# num neighbouring samples to estimate distance upper bound [20]
		if [[ -n ${SC_10X_ARKS_NUMNEIGHBOUR} && ${SC_10X_ARKS_NUMNEIGHBOUR} -gt 0 ]]
		then
			ARKS_OPT="${ARKS_OPT} -B ${SC_10X_ARKS_NUMNEIGHBOUR}"
		else 	
			ARKS_OPT="${ARKS_OPT} -B 20"
		fi
		# -s=FILE output TSV of intra-contig distance/barcode data [disabled]
		if [[ -n ${SC_10X_ARKS_INTRACONTIGTSV} ]]
		then
			ARKS_OPT="${ARKS_OPT} -s ${SC_10X_ARKS_INTRACONTIGTSV}"
		fi
		#-S=FILE output TSV of inter-contig distance/barcode data [disabled]
		if [[ -n ${SC_10X_ARKS_INTERCONTIGTSV} ]]
		then
			ARKS_OPT="${ARKS_OPT} -S ${SC_10X_ARKS_INTERCONTIGTSV}"
		fi
	fi
	# EXTRA OUTPUT OPTIONS	
	if [[ -n ${SC_10X_ARKS_CHECKPOINTS} && ${SC_10X_ARKS_CHECKPOINTS} -ge 0 && ${SC_10X_ARKS_CHECKPOINTS} -le 3 ]]
	then
		ARKS_OPT="${ARKS_OPT} -o ${SC_10X_ARKS_CHECKPOINTS}"
	fi
	# Maximum degree of nodes in graph. All nodes with degree greater than this number will be removed from the graph prior to printing final graph. For no node removal, set to 0 (default: 0)		
	if [[ -n ${SC_10X_ARKS_VERBOSE} && ${SC_10X_ARKS_MAXNODEDEGREE} -ge 0  ]]
	then
		ARKS_OPT="${ARKS_OPT} -d ${SC_10X_ARKS_MAXNODEDEGREE}"
	fi
	
	if [[ -n ${SC_10X_ARKS_THREADS} && ${SC_10X_ARKS_THREADS} -gt 0 ]]
	then
		ARKS_OPT="${ARKS_OPT} -t ${SC_10X_ARKS_THREADS}"
	else
		ARKS_OPT="${ARKS_OPT} -t $(sinfo -p ${SLURM_PARTITION} -o \"%c\" --noheader)"			 
	fi
}

function setLinksOptions()
{
	LINKS_OPT=""
	
	if [[ -n ${SC_10X_LINKS_VERBOSE} && ${SC_10X_LINKS_VERBOSE} -gt 0 ]]
	then
		LINKS_OPT="${LINKS_OPT} -v 1"
	fi
	
	#minimum number of links (k-mer pairs) to compute scaffold (default -l 5, optional)
	if [[ -n ${SC_10X_LINKS_NUMLINKS} && ${SC_10X_LINKS_NUMLINKS} -gt 0 ]]
	then
		LINKS_OPT="${LINKS_OPT} -l ${SC_10X_LINKS_NUMLINKS}"
	else
		LINKS_OPT="${LINKS_OPT} -l 5"
	fi
	
	#maximum link ratio between two best contig pairs (default -a 0.3, optional)
    #     *higher values lead to least accurate scaffolding*
	if [[ -n ${SC_10X_LINKS_MAXLINKRATIO} ]]
	then
		LINKS_OPT="${LINKS_OPT} -a ${SC_10X_LINKS_MAXLINKRATIO}"
	else
		SC_10X_LINKS_MAXLINKRATIO=0.3
		LINKS_OPT="${LINKS_OPT} -a ${SC_10X_LINKS_MAXLINKRATIO}"
	fi
	
	# minimum contig length to consider for scaffolding (default -z 500, optional)
	if [[ -n ${SC_10X_LINKS_MINCONTIGLEN} && ${SC_10X_LINKS_MINCONTIGLEN} -gt 0 ]]
	then
		LINKS_OPT="${LINKS_OPT} -z ${SC_10X_LINKS_MINCONTIGLEN}"
	else
		SC_10X_LINKS_MINCONTIGLEN=10000
		LINKS_OPT="${LINKS_OPT} -z ${SC_10X_LINKS_MINCONTIGLEN}"
	fi
}

#type 0: scaff10x - break10x pipeline		steps: 01_scaff10xprepare, 02_scaff10xbreak10, 03_scaff10xscaff10x, 04_scaff10xbreak10x, 05_scaff10xscaff10x, 06_scaff10xbreak10x, 07_scaff10xStatistics
#type 1: tigmint - arks - links pipeline	steps: 01_arksPrepare, 02_arksLongranger, 03_arksTigmint, 04_arksArks, 05_arksLINKS
#type 2: scaff10x using longranger bam		steps: 01_scaff10xprepare, 02_scaff10xLongrangerAlign, 03_scaff10xPrepareIntermediate, 04_scaff10xScaff10x, 05_scaff10xStatistics
#type 3: break10x using longranger bam		steps: 01_break10xPrepare, 02_break10xLongrangerAlign, 03_break10xPrepareIntermediate, 04_break10xBreak10x, 05_break10xStatistics

myTypes=("01_scaff10xprepare, 02_scaff10xbreak10, 03_scaff10xscaff10x, 04_scaff10xbreak10x, 05_scaff10xscaff10x, 06_scaff10xbreak10x, 07_scaff10xStatistics", 
"01_arksPrepare, 02_arksLongranger, 03_arksTigmint, 04_arksArks", "01_scaff10xprepare, 02_scaff10xLongrangerAlign, 03_scaff10xPrepareIntermediate, 04_scaff10xScaff10x, 05_scaff10xStatistics", 
"01_break10xPrepare, 02_break10xLongrangerAlign, 03_break10xPrepareIntermediate, 04_break10xBreak10x, 05_break10xStatistics") 
if [[ ${SC_10X_TYPE} -eq 0 ]]
then 
    ### 01_scaff10xprepare
    if [[ ${currentStep} -eq 1 ]]
    then
        ### clean up plans 
        for x in $(ls 10x_01_scaff10x*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        if [[ ! -d "${SC_10X_READS}" ]]
        then
        	(>&2 echo "ERROR - set SC_10X_READS to proper 10x read directory")
        	exit 1
   		fi
   		
		numR1Files=0
		for x in ${SC_10X_READS}/${PROJECT_ID}_S*_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz
		do
			if [[ -f ${x} ]]
			then	
				numR1Files=$((${numR1Files}+1))	
			fi
		done
		
		if [[ ${numR1Files} -eq 0 ]]
        then
        	(>&2 echo "ERROR - cannot read 10x R1 files with following pattern: ${SC_10X_READS}/${PROJECT_ID}_S*_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz")
        	exit 1
   		fi
   		
   		numR2Files=0
		for x in ${SC_10X_READS}/${PROJECT_ID}_S*_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz
		do
			if [[ -f ${x} ]]
			then	
				numR2Files=$((${numR2Files}+1))	
			fi
		done
		
		if [[ ${numR2Files} -eq 0 ]]
        then
        	(>&2 echo "ERROR - cannot read 10x R2 files with following pattern: ${SC_10X_READS}/${PROJECT_ID}_S*_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz")
        	exit 1
   		fi
   		
   		if [[ ${numR1Files} -ne ${numR2Files} ]]
        then
        	(>&2 echo "ERROR - 10x R1 files ${numR1Files} does not match R2 files ${numR2Files}")
        	exit 1
   		fi
   		
   		if [[ ! -f ${SC_10X_REF} ]]
   		then
   			(>&2 echo "ERROR - set SC_10X_REF to proper reference fasta file")
        	exit 1	
   		fi
   		
   		echo "if [[ -d ${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID} ]]; then mv ${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID} ${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID}_$(date '+%Y-%m-%d_%H-%M-%S'); fi && mkdir -p ${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID}" > 10x_01_scaff10xprepare_single_${CONT_DB}.${slurmID}.plan
   		echo "mkdir -p ${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID}/ref" >> 10x_01_scaff10xprepare_single_${CONT_DB}.${slurmID}.plan
   		echo "mkdir -p ${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID}/reads" >> 10x_01_scaff10xprepare_single_${CONT_DB}.${slurmID}.plan
   		echo "ln -s -r ${SC_10X_REF} ${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID}/ref" >> 10x_01_scaff10xprepare_single_${CONT_DB}.${slurmID}.plan
   		
   		if [[ -n ${SC_10X_SCAFF10X_READSBC1} && -f ${SC_10X_SCAFF10X_READSBC1} && -n ${SC_10X_SCAFF10X_READSBC2} && -f ${SC_10X_SCAFF10X_READSBC2} ]]
    	then
   			echo "using reads scaff10x_BC_1.fastq scaff10x_BC_2.fastq from a previous run:"
   			echo "${SC_10X_SCAFF10X_READSBC1}"
   			echo "${SC_10X_SCAFF10X_READSBC2}"   			
   		else
	   		for r1 in ${SC_10X_READS}/${PROJECT_ID}_S[0-9]_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz
			do
				id=$(dirname ${r1})
				f1=$(basename ${r1})
				f2=$(echo "${f1}" | sed -e "s:_R1_:_R2_:")
				
				echo "ln -s -f ${id}/${f1} ${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID}/reads"
				echo "ln -s -f ${id}/${f2} ${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID}/reads"										
				echo "echo \"q1=${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID}/reads/${f1}\" >> ${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID}/scaff10x_inputReads.txt"
				echo "echo \"q2=${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID}/reads/${f2}\" >> ${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID}/scaff10x_inputReads.txt"							 
			done >> 10x_01_scaff10xprepare_single_${CONT_DB}.${slurmID}.plan
			
			options="-debug 1 -tmp $(pwd)/${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID}/"
			echo "${SCAFF10X_PATH}/scaff_reads ${options} ${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID}/scaff10x_inputReads.txt scaff10x_BC_1.fastq scaff10x_BC_2.fastq" >> 10x_01_scaff10xprepare_single_${CONT_DB}.${slurmID}.plan
			echo "scaff_reads $(cat ${SCAFF10X_PATH}/version.txt)" > 10x_01_scaff10xprepare_single_${CONT_DB}.${slurmID}.version
		fi
	### 02_scaff10xbreak10x	
	elif [[ ${currentStep} -eq 2 ]]
    then
    	### clean up plans 
        for x in $(ls 10x_02_scaff10*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        setScaff10xOptions
		setBreak10xOptions
    	# add reference 
    	infiles="ref/$(basename ${SC_10X_REF})"
    	# add BCR1 and BCR2 files
    	if [[ -n ${SC_10X_SCAFF10X_READSBC1} && -f ${SC_10X_SCAFF10X_READSBC1} && -n ${SC_10X_SCAFF10X_READSBC2} && -f ${SC_10X_SCAFF10X_READSBC2} ]]
    	then
    		### we need an absolute path if --tmp flag is used in scaff10x 
    		if [[ ! "${SC_10X_SCAFF10X_READSBC1:0:1}" = "/" ]]
    		then 
    			SC_10X_SCAFF10X_READSBC1=$(pwd)/${SC_10X_SCAFF10X_READSBC1}
    		fi
    		
    		if [[ ! "${SC_10X_SCAFF10X_READSBC2:0:1}" = "/" ]]
    		then 
    			SC_10X_SCAFF10X_READSBC2=$(pwd)/${SC_10X_SCAFF10X_READSBC2}
    		fi
    	else
    		SC_10X_SCAFF10X_READSBC1=scaff10x_BC_1.fastq
    		SC_10X_SCAFF10X_READSBC2=scaff10x_BC_2.fastq
    	fi
    	
    	prevExt=$(basename ${SC_10X_REF%.fasta} | awk -F '[_.]' '{print $(NF-1)}')
                
        options="-debug 1 -tmp $(pwd)/${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID}/"
        echo "${SCAFF10X_PATH}/break10x${SCAFF10X_BREAK10X_OPT} ${options} ${infiles} ${SC_10X_SCAFF10X_READSBC1} ${SC_10X_SCAFF10X_READSBC2} ${PROJECT_ID}_${SC_10X_OUTDIR}_${prevExt}b.p.fasta ${PROJECT_ID}_${SC_10X_OUTDIR}_${prevExt}b.p.breaks" > 10x_02_scaff10xscaff10x_single_${CONT_DB}.${slurmID}.plan        
		echo "break10x $(cat ${SCAFF10X_PATH}/version.txt)" >> 10x_02_scaff10xscaff10x_single_${CONT_DB}.${slurmID}.version	
    ### 03_scaff10xscaff10x		
	elif [[ ${currentStep} -eq 3 ]]
    then
    	### clean up plans 
        for x in $(ls 10x_03_scaff10x*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        setScaff10xOptions
		setBreak10xOptions
    	# add reference 
    	infiles="ref/$(basename ${SC_10X_REF})"
    	# add BCR1 and BCR2 files
    	if [[ -n ${SC_10X_SCAFF10X_READSBC1} && -f ${SC_10X_SCAFF10X_READSBC1} && -n ${SC_10X_SCAFF10X_READSBC2} && -f ${SC_10X_SCAFF10X_READSBC2} ]]
    	then
    		### we need an absolute path if --tmp flag is used in scaff10x 
    		if [[ ! "${SC_10X_SCAFF10X_READSBC1:0:1}" = "/" ]]
    		then 
    			SC_10X_SCAFF10X_READSBC1=$(pwd)/${SC_10X_SCAFF10X_READSBC1}
    		fi
    		
    		if [[ ! "${SC_10X_SCAFF10X_READSBC2:0:1}" = "/" ]]
    		then 
    			SC_10X_SCAFF10X_READSBC2=$(pwd)/${SC_10X_SCAFF10X_READSBC2}
    		fi
    	else
    		SC_10X_SCAFF10X_READSBC1=scaff10x_BC_1.fastq
    		SC_10X_SCAFF10X_READSBC2=scaff10x_BC_2.fastq
    	fi
    	
    	prevExt=$(basename ${SC_10X_REF%.fasta} | awk -F '[_.]' '{print $(NF-1)}')
                
        options="-debug 1 -tmp $(pwd)/${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID}/"
        echo "${SCAFF10X_PATH}/scaff10x${SCAFF10X_SCAFF10X_OPT} ${options} ${infiles} ${SC_10X_SCAFF10X_READSBC1} ${SC_10X_SCAFF10X_READSBC2} ${PROJECT_ID}_${SC_10X_OUTDIR}_${prevExt}x.p.fasta" > 10x_03_scaff10xscaff10x_block_${CONT_DB}.${slurmID}.plan
        echo "${SCAFF10X_PATH}/scaff10x${SCAFF10X_SCAFF10X_OPT} ${options} ${PROJECT_ID}_${SC_10X_OUTDIR}_${prevExt}b.p.fasta ${SC_10X_SCAFF10X_READSBC1} ${SC_10X_SCAFF10X_READSBC2} ${PROJECT_ID}_${SC_10X_OUTDIR}_${prevExt}bx.p.fasta" >> 10x_03_scaff10xscaff10x_block_${CONT_DB}.${slurmID}.plan
        
		echo "scaff10x $(cat ${SCAFF10X_PATH}/version.txt)" > 10x_03_scaff10xscaff10x_block_${CONT_DB}.${slurmID}.version
	### 04_scaff10xbreak10x		
	elif [[ ${currentStep} -eq 4 ]]
    then
		### clean up plans 
        for x in $(ls 10x_04_scaff10x*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        setBreak10xOptions
        
        prevExt=$(basename ${SC_10X_REF%.fasta} | awk -F '[_.]' '{print $(NF-1)}')
        # add reference
        inputScaffold1="${PROJECT_ID}_${SC_10X_OUTDIR}_${prevExt}x.p.fasta" 
    	inputScaffold2="${PROJECT_ID}_${SC_10X_OUTDIR}_${prevExt}bx.p.fasta"    	
    	
    	# add BCR1 and BCR2 files
    	if [[ -n ${SC_10X_SCAFF10X_READSBC1} && -f ${SC_10X_SCAFF10X_READSBC1} && -n ${SC_10X_SCAFF10X_READSBC2} && -f ${SC_10X_SCAFF10X_READSBC2} ]]
    	then
    		### we need an absolute path if --tmp flag is used in scaff10x 
    		if [[ ! "${SC_10X_SCAFF10X_READSBC1:0:1}" = "/" ]]
    		then 
    			SC_10X_SCAFF10X_READSBC1=$(pwd)/${SC_10X_SCAFF10X_READSBC1}
    		fi
    		
    		if [[ ! "${SC_10X_SCAFF10X_READSBC2:0:1}" = "/" ]]
    		then 
    			SC_10X_SCAFF10X_READSBC2=$(pwd)/${SC_10X_SCAFF10X_READSBC2}
    		fi
    	else
    		SC_10X_SCAFF10X_READSBC1=scaff10x_BC_1.fastq
    		SC_10X_SCAFF10X_READSBC2=scaff10x_BC_2.fastq
    	fi
                
        options="-debug 1 -tmp $(pwd)/${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID}/"                
        echo "${SCAFF10X_PATH}/break10x${SCAFF10X_BREAK10X_OPT} ${options} ${inputScaffold1} ${SC_10X_SCAFF10X_READSBC1} ${SC_10X_SCAFF10X_READSBC2} ${inputScaffold1%.p.fasta}b.p.fasta ${inputScaffold1%.p.fasta}b.p.breaks" > 10x_04_scaff10xbreak10x_block_${CONT_DB}.${slurmID}.plan
    	echo "${SCAFF10X_PATH}/break10x${SCAFF10X_BREAK10X_OPT} ${options} ${inputScaffold2} ${SC_10X_SCAFF10X_READSBC1} ${SC_10X_SCAFF10X_READSBC2} ${inputScaffold2%.p.fasta}b.p.fasta ${inputScaffold2%.p.fasta}b.p.breaks" >> 10x_04_scaff10xbreak10x_block_${CONT_DB}.${slurmID}.plan      	
        
		echo "break10x $(cat ${SCAFF10X_PATH}/version.txt)" > 10x_04_scaff10xbreak10x_block_${CONT_DB}.${slurmID}.version		
	### 05_scaff10xscaff10x		
	elif [[ ${currentStep} -eq 5 ]]
    then
		### clean up plans 
        for x in $(ls 10x_05_scaff10x*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        setScaff10xOptions
        
        prevExt=$(basename ${SC_10X_REF%.fasta} | awk -F '[_.]' '{print $(NF-1)}')
        # add reference
        inputScaffold1="${PROJECT_ID}_${SC_10X_OUTDIR}_${prevExt}x.p.fasta" 
    	inputScaffold2="${PROJECT_ID}_${SC_10X_OUTDIR}_${prevExt}bx.p.fasta"    	
    	
    	# add BCR1 and BCR2 files
    	if [[ -n ${SC_10X_SCAFF10X_READSBC1} && -f ${SC_10X_SCAFF10X_READSBC1} && -n ${SC_10X_SCAFF10X_READSBC2} && -f ${SC_10X_SCAFF10X_READSBC2} ]]
    	then
    		### we need an absolute path if --tmp flag is used in scaff10x 
    		if [[ ! "${SC_10X_SCAFF10X_READSBC1:0:1}" = "/" ]]
    		then 
    			SC_10X_SCAFF10X_READSBC1=$(pwd)/${SC_10X_SCAFF10X_READSBC1}
    		fi
    		
    		if [[ ! "${SC_10X_SCAFF10X_READSBC2:0:1}" = "/" ]]
    		then 
    			SC_10X_SCAFF10X_READSBC2=$(pwd)/${SC_10X_SCAFF10X_READSBC2}
    		fi
    	else
    		SC_10X_SCAFF10X_READSBC1=scaff10x_BC_1.fastq
    		SC_10X_SCAFF10X_READSBC2=scaff10x_BC_2.fastq
    	fi
          
        options="-debug 1 -tmp $(pwd)/${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID}/"      
        echo "${SCAFF10X_PATH}/scaff10x${SCAFF10X_SCAFF10X_OPT} ${options} ${inputScaffold1} ${SC_10X_SCAFF10X_READSBC1} ${SC_10X_SCAFF10X_READSBC2} ${inputScaffold1%.p.fasta}x.p.fasta" > 10x_05_scaff10xscaff10x_block_${CONT_DB}.${slurmID}.plan
        echo "${SCAFF10X_PATH}/scaff10x${SCAFF10X_SCAFF10X_OPT} ${options} ${inputScaffold1%.p.fasta}b.p.fasta ${SC_10X_SCAFF10X_READSBC1} ${SC_10X_SCAFF10X_READSBC2} ${inputScaffold1%.p.fasta}bx.p.fasta" >> 10x_05_scaff10xscaff10x_block_${CONT_DB}.${slurmID}.plan
        echo "${SCAFF10X_PATH}/scaff10x${SCAFF10X_SCAFF10X_OPT} ${options} ${inputScaffold2%.p.fasta}b.p.fasta ${SC_10X_SCAFF10X_READSBC1} ${SC_10X_SCAFF10X_READSBC2} ${inputScaffold2%.p.fasta}bx.p.fasta" >> 10x_05_scaff10xscaff10x_block_${CONT_DB}.${slurmID}.plan
        
		echo "scaff10x $(cat ${SCAFF10X_PATH}/version.txt)" > 10x_05_scaff10xscaff10x_block_${CONT_DB}.${slurmID}.version
	### 06_scaff10xbreak10x	
	elif [[ ${currentStep} -eq 6 ]]
    then
		### clean up plans 
        for x in $(ls 10x_06_scaff10x*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        setBreak10xOptions
        setScaff10xOptions
        
        prevExt=$(basename ${SC_10X_REF%.fasta} | awk -F '[_.]' '{print $(NF-1)}')
        # add reference
        inputScaffold1="${PROJECT_ID}_${SC_10X_OUTDIR}_${prevExt}x.p.fasta" 
    	inputScaffold2="${PROJECT_ID}_${SC_10X_OUTDIR}_${prevExt}bx.p.fasta"    	
    	
    	# add BCR1 and BCR2 files
    	if [[ -n ${SC_10X_SCAFF10X_READSBC1} && -f ${SC_10X_SCAFF10X_READSBC1} && -n ${SC_10X_SCAFF10X_READSBC2} && -f ${SC_10X_SCAFF10X_READSBC2} ]]
    	then
    		### we need an absolute path if --tmp flag is used in scaff10x 
    		if [[ ! "${SC_10X_SCAFF10X_READSBC1:0:1}" = "/" ]]
    		then 
    			SC_10X_SCAFF10X_READSBC1=$(pwd)/${SC_10X_SCAFF10X_READSBC1}
    		fi
    		
    		if [[ ! "${SC_10X_SCAFF10X_READSBC2:0:1}" = "/" ]]
    		then 
    			SC_10X_SCAFF10X_READSBC2=$(pwd)/${SC_10X_SCAFF10X_READSBC2}
    		fi
    	else
    		SC_10X_SCAFF10X_READSBC1=scaff10x_BC_1.fastq
    		SC_10X_SCAFF10X_READSBC2=scaff10x_BC_2.fastq
    	fi
                
        options="-debug 1 -tmp $(pwd)/${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID}/"                
        echo "${SCAFF10X_PATH}/break10x${SCAFF10X_BREAK10X_OPT} ${options} ${inputScaffold1%.p.fasta}x.p.fasta ${SC_10X_SCAFF10X_READSBC1} ${SC_10X_SCAFF10X_READSBC2} ${inputScaffold1%.p.fasta}xb.p.fasta ${inputScaffold2%.p.fasta}b.p.breaks" > 10x_06_scaff10xbreak10x_block_${CONT_DB}.${slurmID}.plan
        echo "${SCAFF10X_PATH}/break10x${SCAFF10X_BREAK10X_OPT} ${options} ${inputScaffold1%.p.fasta}bx.p.fasta ${SC_10X_SCAFF10X_READSBC1} ${SC_10X_SCAFF10X_READSBC2} ${inputScaffold1%.p.fasta}bxb.p.fasta ${inputScaffold1%.p.fasta}bxb.p.breaks" >> 10x_06_scaff10xbreak10x_block_${CONT_DB}.${slurmID}.plan
        echo "${SCAFF10X_PATH}/break10x${SCAFF10X_BREAK10X_OPT} ${options} ${inputScaffold2%.p.fasta}bx.p.fasta ${SC_10X_SCAFF10X_READSBC1} ${SC_10X_SCAFF10X_READSBC2} ${inputScaffold2%.p.fasta}bxb.p.fasta ${inputScaffold2%.p.fasta}bxb.p.breaks" >> 10x_06_scaff10xbreak10x_block_${CONT_DB}.${slurmID}.plan
        
		echo "break10x $(cat ${SCAFF10X_PATH}/version.txt)" > 10x_06_scaff10xbreak10x_block_${CONT_DB}.${slurmID}.version
	### 07_scaff10xStatistics		
	elif [[ ${currentStep} -eq 7 ]]
    then
		### clean up plans 
        for x in $(ls 10x_07_scaff10x*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
		### run slurm stats - on the master node !!! Because sacct is not available on compute nodes
    	if [[ $(hostname) == "falcon1" || $(hostname) == "falcon2" ]]
        then 
        	bash ${SUBMIT_SCRIPTS_PATH}/slurmStats.sh ${configFile}
    	else
        	cwd=$(pwd)
        	ssh falcon "cd ${cwd} && bash ${SUBMIT_SCRIPTS_PATH}/slurmStats.sh ${configFile}"
    	fi
    	### create assemblyStats plan 
    	echo "${SUBMIT_SCRIPTS_PATH}/assemblyStats.sh ${configFile} 12" > 10x_07_scaff10xStatistics_single_${CONT_DB}.${slurmID}.plan
    	echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > 10x_07_scaff10xStatistics_single_${CONT_DB}.${slurmID}.version	
		echo "$(quast.py --version)" >> 10x_07_scaff10xStatistics_single_${CONT_DB}.${slurmID}.version
    else
        (>&2 echo "step ${currentStep} in SC_10X_TYPE ${SC_10X_TYPE} not supported")
        (>&2 echo "valid steps are: ${myTypes[${SC_10X_TYPE}]}")
        exit 1            
    fi     
elif [[ ${SC_10X_TYPE} -eq 1 ]]
then 
    ### 01_arksPrepare
	### - create directorries, link files
    if [[ ${currentStep} -eq 1 ]]
    then
        ### clean up plans 
        for x in $(ls 10x_01_arks*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        if [[ ! -d "${SC_10X_READS}" ]]
        then
        	(>&2 echo "ERROR - set SC_10X_READS to proper 10x read directory")
        	exit 1
   		fi
   		
		numR1Files=0
		for x in ${SC_10X_READS}/${PROJECT_ID}_S*_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz
		do
			if [[ -f ${x} ]]
			then	
				numR1Files=$((${numR1Files}+1))	
			fi
		done
		
		if [[ ${numR1Files} -eq 0 ]]
        then
        	(>&2 echo "ERROR - cannot read 10x R1 files with following pattern: ${SC_10X_READS}/${PROJECT_ID}_S*_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz")
        	exit 1
   		fi
   		
   		numR2Files=0
		for x in ${SC_10X_READS}/${PROJECT_ID}_S*_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz
		do
			if [[ -f ${x} ]]
			then	
				numR2Files=$((${numR2Files}+1))	
			fi
		done
		
		if [[ ${numR2Files} -eq 0 ]]
        then
        	(>&2 echo "ERROR - cannot read 10x R2 files with following pattern: ${SC_10X_READS}/${PROJECT_ID}_S*_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz")
        	exit 1
   		fi
   		
   		if [[ ${numR1Files} -ne ${numR2Files} ]]
        then
        	(>&2 echo "ERROR - 10x R1 files ${numR1Files} does not match R2 files ${numR2Files}")
        	exit 1
   		fi
   		
   		if [[ ! -f ${SC_10X_REF} ]]
   		then
   			(>&2 echo "ERROR - set SC_10X_REF to proper reference fasta file")
        	exit 1	
   		fi	
   		
   		echo "if [[ -d ${SC_10X_OUTDIR}/arks_${SC_10X_RUNID} ]]; then mv ${SC_10X_OUTDIR}/arks_${SC_10X_RUNID} ${SC_10X_OUTDIR}/arks_${SC_10X_RUNID}_$(date '+%Y-%m-%d_%H-%M-%S'); fi && mkdir -p ${SC_10X_OUTDIR}/arks_${SC_10X_RUNID}" > 10x_01_arksPrepare_single_${CONT_DB}.${slurmID}.plan
   		echo "mkdir -p ${SC_10X_OUTDIR}/arks_${SC_10X_RUNID}/ref" >> 10x_01_arksPrepare_single_${CONT_DB}.${slurmID}.plan
   		echo "mkdir -p ${SC_10X_OUTDIR}/arks_${SC_10X_RUNID}/reads" >> 10x_01_arksPrepare_single_${CONT_DB}.${slurmID}.plan
   		echo "mkdir -p ${SC_10X_OUTDIR}/arks_${SC_10X_RUNID}/tigmint" >> 10x_01_arksPrepare_single_${CONT_DB}.${slurmID}.plan
   		echo "mkdir -p ${SC_10X_OUTDIR}/arks_${SC_10X_RUNID}/arks" >> 10x_01_arksPrepare_single_${CONT_DB}.${slurmID}.plan
   		echo "mkdir -p ${SC_10X_OUTDIR}/arks_${SC_10X_RUNID}/tigmint_arks" >> 10x_01_arksPrepare_single_${CONT_DB}.${slurmID}.plan
   		
   		echo "ln -s -r ${SC_10X_REF} ${SC_10X_OUTDIR}/arks_${SC_10X_RUNID}/ref/${PROJECT_ID}.fasta" >> 10x_01_arksPrepare_single_${CONT_DB}.${slurmID}.plan
   		echo "samtools faidx ${SC_10X_OUTDIR}/arks_${SC_10X_RUNID}/ref/${PROJECT_ID}.fasta" >> 10x_01_arksPrepare_single_${CONT_DB}.${slurmID}.plan
		echo "${PACBIO_BASE_ENV} && bwa index ${SC_10X_OUTDIR}/arks_${SC_10X_RUNID}/ref/${PROJECT_ID}.fasta && ${PACBIO_BASE_ENV_DEACT}" >> 10x_01_arksPrepare_single_${CONT_DB}.${slurmID}.plan
		echo "perl -ne 'chomp; if(/>/){\$ct+=1; print \">\$ct\n\";}else{print \"\$_\n\";} ' < ${SC_10X_OUTDIR}/arks_${SC_10X_RUNID}/ref/${PROJECT_ID}.fasta > ${SC_10X_OUTDIR}/arks_${SC_10X_RUNID}/ref/${PROJECT_ID}.renamed.fasta" >> 10x_01_arksPrepare_single_${CONT_DB}.${slurmID}.plan		
   		
		for r1 in ${SC_10X_READS}/${PROJECT_ID}_S[0-9]_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz
		do
			id=$(dirname ${r1})
			f1=$(basename ${r1})
			f2=$(echo "${f1}" | sed -e "s:_R1_:_R2_:")
			
			echo "ln -s -f ${id}/${f1} ${SC_10X_OUTDIR}/arks_${SC_10X_RUNID}/reads"
			echo "ln -s -f ${id}/${f2} ${SC_10X_OUTDIR}/arks_${SC_10X_RUNID}/reads"										
		done >> 10x_01_arksPrepare_single_${CONT_DB}.${slurmID}.plan		
		
		echo "bwa $(${PACBIO_BASE_ENV} && bwa 2>&1 | grep Version | awk '{print $2}' && ${PACBIO_BASE_ENV_DEACT})" > 10x_01_arksPrepare_single_${CONT_DB}.${slurmID}.version
   		echo "samtools $(${PACBIO_BASE_ENV} && samtools 2>&1 | grep Version | awk '{print $2}' && ${PACBIO_BASE_ENV_DEACT})" >> 10x_01_arksPrepare_single_${CONT_DB}.${slurmID}.version
		
	### 02_arksLongranger
	elif [[ ${currentStep} -eq 2 ]]
    then
        ### clean up plans 
        for x in $(ls 10x_02_arks*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        echo "cd ${SC_10X_OUTDIR}/arks_${SC_10X_RUNID} && ${LONGRANGER_PATH}/longranger basic --id=${PROJECT_ID} --fastqs=reads/" > 10x_02_arksLongranger_single_${CONT_DB}.${slurmID}.plan
        echo "$(${LONGRANGER_PATH}/longranger basic --version | head -n1 | tr -d \"()\")" > 10x_02_arksLongranger_single_${CONT_DB}.${slurmID}.version
    ### 03_arksTigmint 
	elif [[ ${currentStep} -eq 3 ]]
    then
        ### clean up plans 
        for x in $(ls 10x_03_arks*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        setbwaOptions
        setSamtoolsOptions        
        setTigmintOptions
        
        # Align paired-end reads to the draft genome and sort by BX tag.
        echo "${PACBIO_BASE_ENV} && bwa mem${SCAFFOLD_BWA_OPT} ${SC_10X_OUTDIR}/arks_${SC_10X_RUNID}/ref/${PROJECT_ID}.fasta ${SC_10X_OUTDIR}/arks_${SC_10X_RUNID}/${PROJECT_ID}/outs/barcoded.fastq.gz | samtools sort${SCAFFOLD_SAMTOOLS_OPT} -o ${SC_10X_OUTDIR}/arks_${SC_10X_RUNID}/tigmint/${PROJECT_ID}_reads_sortbx.bam && ${PACBIO_BASE_ENV_DEACT}" > 10x_03_arksTigmint_single_${CONT_DB}.${slurmID}.plan
    	# Create molecule extents BED
    	flp="${SC_10X_OUTDIR}/arks_${SC_10X_RUNID}/tigmint/${PROJECT_ID}.as${SC_10X_TIGMINT_MOLECULE_ALNSCORERATIO}.nm${SC_10X_TIGMINT_MOLECULE_MAXMISMATCH}.molecule.size${SC_10X_TIGMINT_MOLECULE_MINMOLSIZE}"    	
    	echo "${TIGMINT_PATH}/tigmint-molecule${TIGMINT_MOLECULE_OPT} ${SC_10X_OUTDIR}/arks_${SC_10X_RUNID}/tigmint/${PROJECT_ID}_reads_sortbx.bam | sort -k1,1 -k2,2n -k3,3n > ${flp}.bed" >> 10x_03_arksTigmint_single_${CONT_DB}.${slurmID}.plan
        # Create molecule extents TSV
        echo "${TIGMINT_PATH}/tigmint-molecule${TIGMINT_MOLECULE_OPT} --tsv -o ${flp}.tsv ${SC_10X_OUTDIR}/arks_${SC_10X_RUNID}/tigmint/${PROJECT_ID}_reads_sortbx.bam" >> 10x_03_arksTigmint_single_${CONT_DB}.${slurmID}.plan
        # Report summary statistics of a Chromium library
        echo "Rscript -e 'rmarkdown::render(\"${TIGMINT_PATH}/../summary.rmd\", \"html_document\", \"${flp}.summary.html\", params= list(input_tsv=\"${flp}.tsv\", output_tsv=\"${flp}.summary.tsv\"))'" >> 10x_03_arksTigmint_single_${CONT_DB}.${slurmID}.plan
        # Compute statistics on the depth of coverage of a BED file.
    	echo "(printf \"Rname\tDepth\tCount\tRsize\tFraction\n\"; awk '\$2 != \$3' ${flp}.bed | bedtools genomecov -g ${SC_10X_OUTDIR}/arks_${SC_10X_RUNID}/ref/${PROJECT_ID}.fasta.fai -i -) > ${flp}.bed.genomecov.tsv" >> 10x_03_arksTigmint_single_${CONT_DB}.${slurmID}.plan
    	# Calculate depth of coverage statistics from bedtools genomecov.
		echo "mlr --tsvlite then filter '\$Rname == \"genome\" && \$Depth > 0' then step -a rsum -f Fraction then put -q '@Depth_count += \$Count; if (is_null(@p25) && \$Fraction_rsum >= 0.25) { @p25 = \$Depth }; if (is_null(@p50) && \$Fraction_rsum >= 0.50) { @p50 = \$Depth }; if (is_null(@p75) && \$Fraction_rsum >= 0.75) { @p75 = \$Depth } end { emitf @Depth_count, @p25, @p50, @p75 }' then rename p25,Depth_p25,p50,Depth_p50,p75,Depth_p75 then put '\$Depth_IQR = \$Depth_p75 - \$Depth_p25' ${flp}.bed.genomecov.tsv > ${flp}.bed.genomecov.stats.tsv" >> 10x_03_arksTigmint_single_${CONT_DB}.${slurmID}.plan
        # Identify breakpoints
		# Make breakpoints BED file
        
        prevExt=$(basename ${SC_10X_REF%.fasta} | awk -F '[_.]' '{print $(NF-1)}')
		cset=$(basename ${SC_10X_REF%.fasta} | awk -F '[_.]' '{print $(NF)}')
		fext="t" ### tigmint
        tigmintOFile=${SC_10X_OUTDIR}/arks_${SC_10X_RUNID}/tigmint/${PROJECT_ID}_${SC_10X_OUTDIR}_${prevExt}${fext}.${cset}.trim${SC_10X_TIGMINT_CUT_TRIM}.window${SC_10X_TIGMINT_CUT_WINDOW}.span${SC_10X_TIGMINT_CUT_SPAN}.breaktigs.fasta
        echo "${TIGMINT_PATH}/tigmint-cut${TIGMINT_CUT_OPT} -o ${tigmintOFile} ${SC_10X_OUTDIR}/arks_${SC_10X_RUNID}/ref/${PROJECT_ID}.fasta ${SC_10X_OUTDIR}/arks_${SC_10X_RUNID}/tigmint/${PROJECT_ID}.as${SC_10X_TIGMINT_MOLECULE_ALNSCORERATIO}.nm${SC_10X_TIGMINT_MOLECULE_MAXMISMATCH}.molecule.size${SC_10X_TIGMINT_MOLECULE_MINMOLSIZE}.bed" >> 10x_03_arksTigmint_single_${CONT_DB}.${slurmID}.plan
		echo "perl -ne 'chomp; if(/>/){\$ct+=1; print \">\$ct\n\";}else{print \"\$_\n\";} ' < ${tigmintOFile} > ${tigmintOFile%.fasta}.renamed.fasta" >> 10x_03_arksTigmint_single_${CONT_DB}.${slurmID}.plan
        
        echo "$(${TIGMINT_PATH}/tigmint-molecule --version)" > 10x_03_arksTigmint_single_${CONT_DB}.${slurmID}.version
		echo "$(${TIGMINT_PATH}/tigmint-cut --version)" >> 10x_03_arksTigmint_single_${CONT_DB}.${slurmID}.version
		echo "mlr $(mlr --version | awk '{print $2}') ">> 10x_03_arksTigmint_single_${CONT_DB}.${slurmID}.version
	### 04_arksArks
	elif [[ ${currentStep} -eq 4 ]]
    then
        ### clean up plans 
        for x in $(ls 10x_04_arks*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done

		setArksOptions
		echo "echo ${SC_10X_OUTDIR}/arks_${SC_10X_RUNID}/${PROJECT_ID}/outs/barcoded.fastq.gz > ${SC_10X_OUTDIR}/arks_${SC_10X_RUNID}/longrangerReads.fof" > 10x_04_arksArks_single_${CONT_DB}.${slurmID}.plan
    	echo "${ARKS_PATH}/calcBarcodeMultiplicities.pl ${SC_10X_OUTDIR}/arks_${SC_10X_RUNID}/longrangerReads.fof > ${SC_10X_OUTDIR}/arks_${SC_10X_RUNID}/longrangerReads_multiplicities.csv" >> 10x_04_arksArks_single_${CONT_DB}.${slurmID}.plan
    	
    	## run type
		# -p can be one of:
	    #                    1) full    uses the full ARKS process (kmerize draft, kmerize and align chromium reads, scaffold).
	    #                    2) align   skips kmerizing of draft and starts with kmerizing and aligning chromium reads.
	    #                    3) graph   skips kmerizing draft and kmerizing/aligning chromium reads and only scaffolds.
		if [[ -n ${SC_10X_ARKS_RUNTYPE} ]]
		then
			if [[ "x${SC_10X_ARKS_RUNTYPE}" != "xfull" && "x${SC_10X_ARKS_RUNTYPE}" != "xalign" && "x${SC_10X_ARKS_RUNTYPE}" != "xgraph" ]]
			then 
				SC_10X_ARKS_RUNTYPE="full"
			fi	
			ARKS_OPT="${ARKS_OPT} -p ${SC_10X_ARKS_RUNTYPE}"
		fi
		
		prevExt=$(basename ${SC_10X_REF%.fasta} | awk -F '[_.]' '{print $(NF-1)}')
		cset=$(basename ${SC_10X_REF%.fasta} | awk -F '[_.]' '{print $(NF)}')
		fext="t" ### tigmint
		tigmintOFile=${SC_10X_OUTDIR}/arks_${SC_10X_RUNID}/tigmint/${PROJECT_ID}_${SC_10X_OUTDIR}_${prevExt}${fext}.${cset}.trim${SC_10X_TIGMINT_CUT_TRIM}.window${SC_10X_TIGMINT_CUT_WINDOW}.span${SC_10X_TIGMINT_CUT_SPAN}.breaktigs.renamed.fasta
        bopt=${SC_10X_OUTDIR}/arks_${SC_10X_RUNID}/tigmint_arks/${PROJECT_ID}_${SC_10X_OUTDIR}_${prevExt}${fext}.${cset}_c${SC_10X_ARKS_MINREADPAIRS}_m${SC_10X_ARKS_MULTIPLICITY}_k${SC_10X_ARKS_KMER}_r${SC_10X_ARKS_PVALUE}_e${SC_10X_ARKS_ENDLEN}_z${SC_10X_ARKS_MINCONTIGLEN}
        
		# ARCS
		# Create a graph of linked contigs using ARCS.
		
		if [[ "${SC_10X_ARKS_RUNTYPE}" == full ]]
		then 
			echo "${ARKS_PATH}/arks ${ARKS_OPT} -a ${SC_10X_OUTDIR}/arks_${SC_10X_RUNID}/longrangerReads_multiplicities.csv -f ${tigmintOFile} -b ${bopt} ${SC_10X_OUTDIR}/arks_${SC_10X_RUNID}/${PROJECT_ID}/outs/barcoded.fastq.gz"
			# Generate TSV from ARKS
			echo "python ${ARKS_PATH}/../Examples/makeTSVfile.py ${bopt} ${bopt%_original.gv}.tigpair_checkpoint.tsv ${tigmintOFile}"		
		else 
			(>&2 echo "SC_10X_ARKS_RUNTYPE: ${SC_10X_ARKS_RUNTYPE} not supported yet")
	    	exit 1
		fi >> 10x_04_arksArks_single_${CONT_DB}.${slurmID}.plan   
		
		 	    	
    	echo "$(${ARKS_PATH}/arks --version | grep VERSION | awk '{print $3}')" > 10x_04_arksArks_single_${CONT_DB}.${slurmID}.version
   	### 05_arksLinks
	elif [[ ${currentStep} -eq 5 ]]
    then
        ### clean up plans 
        for x in $(ls 10x_05_arks*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
    
    	if [[ -f ${SC_10X_OUTDIR}/arks_${SC_10X_RUNID}/tigmint_arks/empty.fof ]]; 
    	then 
    		echo "mv ${SC_10X_OUTDIR}/arks_${SC_10X_RUNID}/tigmint_arks/empty.fof ${SC_10X_OUTDIR}/arks_${SC_10X_RUNID}/tigmint_arks/empty_$(date '+%Y-%m-%d_%H-%M-%S').fof"; 
    		echo "touch ${SC_10X_OUTDIR}/arks_${SC_10X_RUNID}/tigmint_arks/empty.fof"
    	else
    		echo "touch ${SC_10X_OUTDIR}/arks_${SC_10X_RUNID}/tigmint_arks/empty.fof"
		fi > 10x_05_arksLinks_single_${CONT_DB}.${slurmID}.plan
		
		setLinksOptions
		
		prevExt=$(basename ${SC_10X_REF%.fasta} | awk -F '[_.]' '{print $(NF-1)}')
		cset=$(basename ${SC_10X_REF%.fasta} | awk -F '[_.]' '{print $(NF)}')
		fext="t" ### tigmint
		tigmintOFile=${SC_10X_OUTDIR}/arks_${SC_10X_RUNID}/tigmint/${PROJECT_ID}_${SC_10X_OUTDIR}_${prevExt}${fext}.${cset}.trim${SC_10X_TIGMINT_CUT_TRIM}.window${SC_10X_TIGMINT_CUT_WINDOW}.span${SC_10X_TIGMINT_CUT_SPAN}.breaktigs.renamed.fasta
		
		empt=${SC_10X_OUTDIR}/arks_${SC_10X_RUNID}/tigmint_arks/empty.fof
		scaf=${SC_10X_OUTDIR}/arks_${SC_10X_RUNID}/tigmint_arks/${PROJECT_ID}_${SC_10X_OUTDIR}_${prevExt}${fext}.${cset}.trim${SC_10X_TIGMINT_CUT_TRIM}.window${SC_10X_TIGMINT_CUT_WINDOW}.span${SC_10X_TIGMINT_CUT_SPAN}_z${SC_10X_LINKS_MINCONTIGLEN}_l${SC_10X_LINKS_MINCONTIGLEN}_a${SC_10X_LINKS_MAXLINKRATIO}.scaffolds.fa
		echo "${LINKS_PATH}/LINKS -f ${tigmintOFile} -s ${empt} -b ${scaf} ${LINKS_OPT}" >> 10x_05_arksLinks_single_${CONT_DB}.${slurmID}.plan
   	else
        (>&2 echo "step ${currentStep} in SC_10X_TYPE ${SC_10X_TYPE} not supported")
        (>&2 echo "valid steps are: ${myTypes[${SC_10X_TYPE}]}")
        exit 1            
    fi
#type 2: scaff10x using longranger bam		steps: 01_scaff10xprepare, 02_scaff10xLongrangerAlign, 03_scaff10xPrepareIntermediate, 04_scaff10xScaff10x, 05_scaff10xStatistics    
elif [[ ${SC_10X_TYPE} -eq 2 ]]
then 
    ### 01_scaff10xprepare
    if [[ ${currentStep} -eq 1 ]]
    then
        ### clean up plans 
        for x in $(ls 10x_01_scaff10x*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        if [[ ! -d "${SC_10X_READS}" ]]
        then
        	(>&2 echo "ERROR - set SC_10X_READS to proper 10x read directory")
        	exit 1
   		fi
   		
   		if [[ ! -f ${SC_10X_REF} ]]
   		then
   			(>&2 echo "ERROR - set SC_10X_REF to proper reference fasta file")
        	exit 1	
   		fi
   		   		   		
   		echo "if [[ -d ${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID} ]]; then mv ${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID} ${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID}_$(date '+%Y-%m-%d_%H-%M-%S'); fi && mkdir -p ${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID}" > 10x_01_scaff10xprepare_single_${CONT_DB}.${slurmID}.plan
   		echo "mkdir -p ${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID}/ref" >> 10x_01_scaff10xprepare_single_${CONT_DB}.${slurmID}.plan
   		echo "mkdir -p ${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID}/bams" >> 10x_01_scaff10xprepare_single_${CONT_DB}.${slurmID}.plan
   		  		
   		REFNAME=$(basename ${SC_10X_REF})
   		if [[ -f ${SC_10X_REF_EXCLUDELIST} ]]
   		then
   			echo "${SEQKIT_PATH} grep -v -f ${SC_10X_REF_EXCLUDELIST} ${SC_10X_REF} > ${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID}/ref/${REFNAME}"
   			echo "${SEQKIT_PATH} grep -f ${SC_10X_REF_EXCLUDELIST} ${SC_10X_REF} > ${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID}/ref/exclude.fasta"
   		else
   			echo "ln -s -r -f ${SC_10X_REF} ${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID}/ref/${REFNAME}" 
   		fi >> 10x_01_scaff10xprepare_single_${CONT_DB}.${slurmID}.plan   		
   		
   		## we need a special reference fastq file with other names  
   		echo "${SCAFF10X_PATH}/scaff-bin/scaff_fastq -name tarseq -len 10 ${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID}/ref/${REFNAME} ${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID}/ref/tarseq.fastq ${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID}/ref/tarseq.tag" >> 10x_01_scaff10xprepare_single_${CONT_DB}.${slurmID}.plan
   		## convert the reference fastq into fasta 
   		echo "sed -n '1~4s/^@/>/p;2~4p' ${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID}/ref/tarseq.fastq > ${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID}/ref/tarseq.fasta" >> 10x_01_scaff10xprepare_single_${CONT_DB}.${slurmID}.plan

		echo "scaff_fastq $(${SCAFF10X_PATH}/scaff10x | grep Version)" > 10x_01_scaff10xprepare_single_${CONT_DB}.${slurmID}.version
	## 02_scaff10xLongrangerAlign
    elif [[ ${currentStep} -eq 2 ]]
    then
        ### clean up plans 
        for x in $(ls 10x_02_scaff10x*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
                
        if [[ ! -f "${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID}/ref/tarseq.fasta" ]]
        then
    		(>&2 echo "ERROR - cannot find reference fasta file ${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID}/ref/tarseq.fasta. Rerun step 01_scaff10xprepare!")
        	exit 1
   		fi
        
		if [[ ! -d ${TENX_PATH} ]]
        then 
        	(>&2 echo "ERROR - cannot find 10x reads. Variable TENX_PATH has to be set to a directoty containing 10x reads.")
        	exit 1
    	fi
    	 
        if [[ ! -d ${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID}/ref/refdata-tarseq ]]
        then
        	echo "cd ${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID}/ref && ${LONGRANGER_PATH}/longranger mkref tarseq.fasta && cd ../../../ " 
        	echo "cd ${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID}/bams && ${LONGRANGER_PATH}/longranger align --id=10x_${PROJECT_ID}_longrangerAlign --fastqs=${TENX_PATH} --sample=${PROJECT_ID} --reference=../ref/refdata-tarseq --jobmode=slurm --localcores=24 --localmem=128 --maxjobs=500 --jobinterval=5000 --disable-ui --nopreflight && cd ../../../"
    	else 
    		(>&2 echo "[WARNING] Using previously created reference file ${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID}/ref/refdata-tarseq. Please remove that folder to rerun longranger mkref" )
    		echo "cd ${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID}/bams && ${LONGRANGER_PATH}/longranger align --id=10x_${PROJECT_ID}_longrangerAlign --fastqs=${TENX_PATH} --sample=${PROJECT_ID} --reference=../ref/refdata-tarseq --jobmode=slurm --localcores=24 --localmem=128 --maxjobs=500 --jobinterval=5000 --disable-ui --nopreflight && cd ../../../"
    	fi > 10x_02_scaff10xLongrangerAlign_single_${CONT_DB}.${slurmID}.plan                
        
        echo "$(${LONGRANGER_PATH}/longranger mkref --version | head -n1)" > 10x_02_scaff10xLongrangerAlign_single_${CONT_DB}.${slurmID}.version
        echo "$(${LONGRANGER_PATH}/longranger align --version | head -n1)" >> 10x_02_scaff10xLongrangerAlign_single_${CONT_DB}.${slurmID}.version
    
    ## 03_scaff10xPrepareIntermediate
    elif [[ ${currentStep} -eq 3 ]]
    then
        ### clean up plans 
        for x in $(ls 10x_03_scaff10x*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        bam=${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID}/bams/10x_${PROJECT_ID}_longrangerAlign/outs/possorted_bam.bam
   		if [[ ! -f ${bam} ]]
   		then 
   			(>&2 echo "ERROR - cannot final duplicate marked bam file: ${bam}!")
        	exit 1
   		fi 
   		
   		# we have to find the BX tag and append it to the fasta name 
   		# unfortunatelly the column is not fixed, in my example data its 19, 21 and 20 --> to be conservative: loop through all optional fields 12 - NF
   		
   		if [[ -n ${SC_10X_SCAFF10X_LONGREAD} && ${SC_10X_SCAFF10X_LONGREAD} -ne 0 ]]
		then
			echo "samtools view ${bam}  | awk '(\$4!=0)&&(\$2<100)&&(\$5>=0){ for(i = 12; i <= NF; i++) { if (\$i ~ \"BX:Z\") {print \$1\"_\"substr(\$i,6,16),\$2,\$3,\$4,\$5; break; } } }' > ${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID}/bams/dalign.dat"
		else
			echo "samtools view ${bam}  | awk '(\$4!=0)&&(\$2<100)&&(\$5>0){ for(i = 12; i <= NF; i++) { if (\$i ~ \"BX:Z\") {print \$1\"_\"substr(\$i,6,16),\$2,\$3,\$4,\$5; break; } } }' > ${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID}/bams/dalign.dat"
		fi > 10x_03_scaff10xPrepareIntermediate_single_${CONT_DB}.${slurmID}.plan 
   		
	## 04_scaff10xScaff10x   		
   	elif [[ ${currentStep} -eq 4 ]]
    then
        ### clean up plans 
        for x in $(ls 10x_04_scaff10x*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        dat=${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID}/bams/dalign.dat
        if [[ ! -f ${dat} ]]
   		then 
   			(>&2 echo "ERROR - dat file missing : ${dat}!. Rerun step 03_scaff10xPrepareIntermediate!")
        	exit 1
   		fi
   		
   		REFNAME=$(basename ${SC_10X_REF})
   		
   		if [[ ! -f ${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID}/${REFNAME} ]]
   		then
   			ln -s -r ${SC_10X_REF} ${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID}/${REFNAME}
   		fi > 10x_04_scaff10xScaff10x_single_${CONT_DB}.${slurmID}.plan
   		
   		setScaff10xOptions
		
    	prevExt=$(basename ${SC_10X_REF%.fasta} | awk -F '[_.]' '{print $(NF-1)}')
        ext=$(basename ${SC_10X_REF%.fasta} | awk -F '[_.]' '{print $(NF)}')
                
        options="-debug 1 -dat ../bams/dalign.dat"
        echo "cd ${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID}/ && ${SCAFF10X_PATH}/scaff10x${SCAFF10X_SCAFF10X_OPT} ${options} ${REFNAME} ${PROJECT_ID}_${SC_10X_OUTDIR}_${prevExt}x.${ext}.fasta && cd  ../../" >> 10x_04_scaff10xScaff10x_single_${CONT_DB}.${slurmID}.plan        
		echo "scaff10x $(${SCAFF10X_PATH}/scaff10x | grep Version)" > 10x_04_scaff10xScaff10x_single_${CONT_DB}.${slurmID}.version
	## 05_scaff10xStatistics   		
   	elif [[ ${currentStep} -eq 5 ]]
    then
        ### clean up plans 
        for x in $(ls 10x_05_scaff10x*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done	
        
        ### run slurm stats - on the master node !!! Because sacct is not available on compute nodes
    	if [[ $(hostname) == "falcon1" || $(hostname) == "falcon2" ]]
        then 
        	bash ${SUBMIT_SCRIPTS_PATH}/slurmStats.sh ${configFile}
    	else
        	cwd=$(pwd)
        	ssh falcon "cd ${cwd} && bash ${SUBMIT_SCRIPTS_PATH}/slurmStats.sh ${configFile}"
    	fi
    	### create assemblyStats plan 
    	echo "${SUBMIT_SCRIPTS_PATH}/assemblyStats.sh ${configFile} 12" > 10x_05_scaff10xStatistics_single_${CONT_DB}.${slurmID}.plan
    	echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > 10x_05_scaff10xStatistics_single_${CONT_DB}.${slurmID}.version	
		echo "$(quast.py --version)" >> 10x_05_scaff10xStatistics_single_${CONT_DB}.${slurmID}.version	
   	else   		
        (>&2 echo "step ${currentStep} in SC_10X_TYPE ${SC_10X_TYPE} not supported")
        (>&2 echo "valid steps are: ${myTypes[${SC_10X_TYPE}]}")
        exit 1            
    fi	
#type 3: break10x using longranger bam		steps: 01_break10xPrepare, 02_break10xLongrangerAlign, 03_break10xPrepareIntermediate, 04_break10xBreak10x, 05_break10xStatistics    
elif [[ ${SC_10X_TYPE} -eq 3 ]]
then 
    ### 01_break10xPrepare
    if [[ ${currentStep} -eq 1 ]]
    then
        ### clean up plans 
        for x in $(ls 10x_01_break10x*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        if [[ ! -d "${SC_10X_READS}" ]]
        then
        	(>&2 echo "ERROR - set SC_10X_READS to proper 10x read directory")
        	exit 1
   		fi
   		
   		if [[ ! -f ${SC_10X_REF} ]]
   		then
   			(>&2 echo "ERROR - set SC_10X_REF to proper reference fasta file")
        	exit 1	
   		fi
   		   		   		
   		echo "if [[ -d ${SC_10X_OUTDIR}/break10x_${SC_10X_RUNID} ]]; then mv ${SC_10X_OUTDIR}/break10x_${SC_10X_RUNID} ${SC_10X_OUTDIR}/break10x_${SC_10X_RUNID}_$(date '+%Y-%m-%d_%H-%M-%S'); fi && mkdir -p ${SC_10X_OUTDIR}/break10x_${SC_10X_RUNID}" > 10x_01_break10xPrepare_single_${CONT_DB}.${slurmID}.plan
   		echo "mkdir -p ${SC_10X_OUTDIR}/break10x_${SC_10X_RUNID}/ref" >> 10x_01_break10xPrepare_single_${CONT_DB}.${slurmID}.plan
   		echo "mkdir -p ${SC_10X_OUTDIR}/break10x_${SC_10X_RUNID}/bams" >> 10x_01_break10xPrepare_single_${CONT_DB}.${slurmID}.plan
   		
   		REFNAME=$(basename ${SC_10X_REF})
   		if [[ -f ${SC_10X_REF_EXCLUDELIST} ]]
   		then
   			echo "${SEQKIT_PATH} grep -v -f ${SC_10X_REF_EXCLUDELIST} ${SC_10X_REF} > ${SC_10X_OUTDIR}/break10x_${SC_10X_RUNID}/ref/${REFNAME}"
   			echo "${SEQKIT_PATH} grep -f ${SC_10X_REF_EXCLUDELIST} ${SC_10X_REF} > ${SC_10X_OUTDIR}/break10x_${SC_10X_RUNID}/ref/exclude.fasta"
   		else
   			echo "ln -s -r -f  ${SC_10X_REF} ${SC_10X_OUTDIR}/break10x_${SC_10X_RUNID}/ref/${REFNAME}" 
   		fi >> 10x_01_break10xPrepare_single_${CONT_DB}.${slurmID}.plan   		
   		
   		## we need a special reference fastq file with other names  
   		echo "${SCAFF10X_PATH}/scaff-bin/scaff_fastq -name tarseq -len 10 ${SC_10X_OUTDIR}/break10x_${SC_10X_RUNID}/ref/${REFNAME} ${SC_10X_OUTDIR}/break10x_${SC_10X_RUNID}/ref/tarseq.fastq ${SC_10X_OUTDIR}/break10x_${SC_10X_RUNID}/ref/tarseq.tag" >> 10x_01_break10xPrepare_single_${CONT_DB}.${slurmID}.plan
   		## convert the reference fastq into fasta 
   		echo "sed -n '1~4s/^@/>/p;2~4p' ${SC_10X_OUTDIR}/break10x_${SC_10X_RUNID}/ref/tarseq.fastq > ${SC_10X_OUTDIR}/break10x_${SC_10X_RUNID}/ref/tarseq.fasta" >> 10x_01_break10xPrepare_single_${CONT_DB}.${slurmID}.plan

		echo "scaff_fastq scaff10x $(${SCAFF10X_PATH}/scaff10x | grep Version)" > 10x_01_break10xPrepare_single_${CONT_DB}.${slurmID}.version
	## 02_break10xLongrangerAlign
    elif [[ ${currentStep} -eq 2 ]]
    then
        ### clean up plans 
        for x in $(ls 10x_02_break10x*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
                
        if [[ ! -f "${SC_10X_OUTDIR}/break10x_${SC_10X_RUNID}/ref/tarseq.fasta" ]]
        then
    		(>&2 echo "ERROR - cannot find reference fasta file ${SC_10X_OUTDIR}/break10x_${SC_10X_RUNID}/ref/tarseq.fasta. Rerun step 01_break10xPrepare!")
        	exit 1
   		fi
        
		if [[ ! -d ${TENX_PATH} ]]
        then 
        	(>&2 echo "ERROR - cannot find 10x reads. Variable TENX_PATH has to be set to a directoty containing 10x reads.")
        	exit 1
    	fi
    	 
        if [[ ! -d ${SC_10X_OUTDIR}/break10x_${SC_10X_RUNID}/ref/refdata-tarseq ]]
        then
        	echo "cd ${SC_10X_OUTDIR}/break10x_${SC_10X_RUNID}/ref && ${LONGRANGER_PATH}/longranger mkref tarseq.fasta && cd ../../../ " 
        	echo "cd ${SC_10X_OUTDIR}/break10x_${SC_10X_RUNID}/bams && ${LONGRANGER_PATH}/longranger align --id=10x_${PROJECT_ID}_longrangerAlign --fastqs=${TENX_PATH} --sample=${PROJECT_ID} --reference=../ref/refdata-tarseq --jobmode=slurm --localcores=24 --localmem=128 --maxjobs=500 --jobinterval=5000 --disable-ui --nopreflight && cd ../../../"
    	else 
    		(>&2 echo "[WARNING] Using previously created reference file ${SC_10X_OUTDIR}/break10x_${SC_10X_RUNID}/ref/refdata-tarseq. Please remove that folder to rerun longranger mkref from scratch" )
    		echo "cd ${SC_10X_OUTDIR}/break10x_${SC_10X_RUNID}/bams && ${LONGRANGER_PATH}/longranger align --id=10x_${PROJECT_ID}_longrangerAlign --fastqs=${TENX_PATH} --sample=${PROJECT_ID} --reference=../ref/refdata-tarseq --jobmode=slurm --localcores=24 --localmem=128 --maxjobs=500 --jobinterval=5000 --disable-ui --nopreflight && cd ../../../"
		fi > 10x_02_break10xLongrangerAlign_single_${CONT_DB}.${slurmID}.plan                
        
        echo "$(${LONGRANGER_PATH}/longranger mkref --version | head -n1)" > 10x_02_break10xLongrangerAlign_single_${CONT_DB}.${slurmID}.version
        echo "$(${LONGRANGER_PATH}/longranger align --version | head -n1)" >> 10x_02_break10xLongrangerAlign_single_${CONT_DB}.${slurmID}.version
    
    ## 03_break10xPrepareIntermediate
    elif [[ ${currentStep} -eq 3 ]]
    then
        ### clean up plans 
        for x in $(ls 10x_03_break10x*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        bam=${SC_10X_OUTDIR}/break10x_${SC_10X_RUNID}/bams/10x_${PROJECT_ID}_longrangerAlign/outs/possorted_bam.bam
   		if [[ ! -f ${bam} ]]
   		then 
   			(>&2 echo "ERROR - cannot final duplicate marked bam file: ${bam}!")
        	exit 1
   		fi 
   		
   		# we have to find the BX tag and append it to the fasta name 
   		# unfortunatelly the column is not fixed, in my example data its 19, 21 and 20 --> to be conservative: loop through all optional fields 12 - NF
   		
   		if [[ -n ${SC_10X_SCAFF10X_LONGREAD} && ${SC_10X_SCAFF10X_LONGREAD} -ne 0 ]]
		then
			echo "samtools view ${bam}  | awk '(\$4!=0)&&(\$2<100)&&(\$5>=0){ for(i = 12; i <= NF; i++) { if (\$i ~ \"BX:Z\") {print \$1\"_\"substr(\$i,6,16),\$2,\$3,\$4,\$5; break; } } }' > ${SC_10X_OUTDIR}/break10x_${SC_10X_RUNID}/bams/dalign.dat"
		else
			echo "samtools view ${bam}  | awk '(\$4!=0)&&(\$2<100)&&(\$5>0){ for(i = 12; i <= NF; i++) { if (\$i ~ \"BX:Z\") {print \$1\"_\"substr(\$i,6,16),\$2,\$3,\$4,\$5; break; } } }' > ${SC_10X_OUTDIR}/break10x_${SC_10X_RUNID}/bams/dalign.dat"
		fi > 10x_03_break10xPrepareIntermediate_single_${CONT_DB}.${slurmID}.plan        

	## 04_break10xScaff10x   		
   	elif [[ ${currentStep} -eq 4 ]]
    then
        ### clean up plans 
        for x in $(ls 10x_04_break10x*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        dat=${SC_10X_OUTDIR}/break10x_${SC_10X_RUNID}/bams/dalign.dat
        if [[ ! -f ${dat} ]]
   		then 
   			(>&2 echo "ERROR - dat file missing : ${dat}!. Rerun step 03_break10xPrepareIntermediate!")
        	exit 1
   		fi
   		
   		REFNAME=$(basename ${SC_10X_REF})
   		
   		if [[ ! -f ${SC_10X_OUTDIR}/break10x_${SC_10X_RUNID}/${REFNAME} ]]
   		then
   			ln -s -r ${SC_10X_REF} ${SC_10X_OUTDIR}/break10x_${SC_10X_RUNID}/${REFNAME}
		fi > 10x_04_break10xBreak10x_single_${CONT_DB}.${slurmID}.plan
   		
   		setBreak10xOptions
		
    	prevExt=$(basename ${SC_10X_REF%.fasta} | awk -F '[_.]' '{print $(NF-1)}')
        ext=$(basename ${SC_10X_REF%.fasta} | awk -F '[_.]' '{print $(NF)}')
                
        options="-debug 1 -dat ../bams/dalign.dat"
        echo "cd ${SC_10X_OUTDIR}/break10x_${SC_10X_RUNID}/ && ${SCAFF10X_PATH}/break10x${SCAFF10X_BREAK10X_OPT} ${options} ${REFNAME} ${PROJECT_ID}_${SC_10X_OUTDIR}_${prevExt}b.${ext}.fasta ${PROJECT_ID}_${SC_10X_OUTDIR}_${prevExt}b.${ext}.breaks && cd  ../../" >> 10x_04_break10xBreak10x_single_${CONT_DB}.${slurmID}.plan        
		echo "break10x $(${SCAFF10X_PATH}/break10x | grep Version)" > 10x_04_break10xBreak10x_single_${CONT_DB}.${slurmID}.version
	## 05_break10xStatistics   		
   	elif [[ ${currentStep} -eq 5 ]]
    then
        ### clean up plans 
        for x in $(ls 10x_05_break10x*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done	
        
        ### run slurm stats - on the master node !!! Because sacct is not available on compute nodes
    	if [[ $(hostname) == "falcon1" || $(hostname) == "falcon2" ]]
        then 
        	bash ${SUBMIT_SCRIPTS_PATH}/slurmStats.sh ${configFile}
    	else
        	cwd=$(pwd)
        	ssh falcon "cd ${cwd} && bash ${SUBMIT_SCRIPTS_PATH}/slurmStats.sh ${configFile}"
    	fi
    	### create assemblyStats plan 
    	echo "${SUBMIT_SCRIPTS_PATH}/assemblyStats.sh ${configFile} 12" > 10x_05_break10xStatistics_single_${CONT_DB}.${slurmID}.plan
    	echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > 10x_05_break10xStatistics_single_${CONT_DB}.${slurmID}.version	
		echo "$(quast.py --version)" >> 10x_05_break10xStatistics_single_${CONT_DB}.${slurmID}.version                
   	else
        (>&2 echo "step ${currentStep} in SC_10X_TYPE ${SC_10X_TYPE} not supported")
        (>&2 echo "valid steps are: ${myTypes[${SC_10X_TYPE}]}")
        exit 1            
    fi    
else
    (>&2 echo "unknown SC_10X_TYPE ${SC_10X_TYPE}")
    (>&2 echo "supported types")
    x=0; while [ $x -lt ${#myTypes[*]} ]; do (>&2 echo "${myTypes[${x}]}"); done 
    exit 1
fi

exit 0