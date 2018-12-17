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

if [[ -z ${SCAFF10X_PATH} || ! -f ${SCAFF10X_PATH}/scaff_reads ]]
then
	(>&2 echo "Variable SCAFF10X_PATH must be set to a proper scaff10x installation directory!!")
    exit 1
fi

function setScaff10xOptions()
{
	SCAFF10X_SCAFF10X_OPT=""
	
	if [[ -n ${SCAF_SCAFF10X_SCAFF10X_THREADS} && ${SCAF_SCAFF10X_SCAFF10X_THREADS} -ne 0 ]]
	then
		SCAFF10X_SCAFF10X_OPT="${SCAFF10X_SCAFF10X_OPT} -nodes ${SCAF_SCAFF10X_SCAFF10X_THREADS}"		
	fi
	
	if [[ -n ${SCAF_SCAFF10X_SCAFF10X_ALIGNER} ]]
	then
		SCAFF10X_SCAFF10X_OPT="${SCAFF10X_SCAFF10X_OPT} -align ${SCAF_SCAFF10X_SCAFF10X_ALIGNER}"		
	fi
	
	if [[ -n ${SCAF_SCAFF10X_SCAFF10X_SCORE} && ${SCAF_SCAFF10X_SCAFF10X_SCORE} -ne 0 ]]
	then
		SCAFF10X_SCAFF10X_OPT="${SCAFF10X_SCAFF10X_OPT} -score ${SCAF_SCAFF10X_SCAFF10X_SCORE}"		
	fi
	
	if [[ -n ${SCAF_SCAFF10X_SCAFF10X_MATRIX} && ${SCAF_SCAFF10X_SCAFF10X_MATRIX} -ne 0 ]]
	then
		SCAFF10X_SCAFF10X_OPT="${SCAFF10X_SCAFF10X_OPT} -matrix ${SCAF_SCAFF10X_SCAFF10X_MATRIX}"		
	fi
	
	if [[ -n ${SCAF_SCAFF10X_SCAFF10X_MINREADS} && ${SCAF_SCAFF10X_SCAFF10X_MINREADS} -ne 0 ]]
	then
		SCAFF10X_SCAFF10X_OPT="${SCAFF10X_SCAFF10X_OPT} -reads ${SCAF_SCAFF10X_SCAFF10X_MINREADS}"		
	fi
	
	if [[ -n ${SCAF_SCAFF10X_SCAFF10X_LONGREAD} && ${SCAF_SCAFF10X_SCAFF10X_LONGREAD} -ne 0 ]]
	then
		SCAFF10X_SCAFF10X_OPT="${SCAFF10X_SCAFF10X_OPT} -longread ${SCAF_SCAFF10X_SCAFF10X_LONGREAD}"		
	fi
	
	if [[ -n ${SCAF_SCAFF10X_SCAFF10X_GAPSIZE} && ${SCAF_SCAFF10X_SCAFF10X_GAPSIZE} -ne 0 ]]
	then
		SCAFF10X_SCAFF10X_OPT="${SCAFF10X_SCAFF10X_OPT} -gap ${SCAF_SCAFF10X_SCAFF10X_GAPSIZE}"		
	fi
	
	if [[ -n ${SCAF_SCAFF10X_SCAFF10X_EDGELEN} && ${SCAF_SCAFF10X_SCAFF10X_EDGELEN} -ne 0 ]]
	then
		SCAFF10X_SCAFF10X_OPT="${SCAFF10X_SCAFF10X_OPT} -edge ${SCAF_SCAFF10X_SCAFF10X_EDGELEN}"		
	fi
	
	if [[ -n ${SCAF_SCAFF10X_SCAFF10X_MINSHAREDBARCODES} && ${SCAF_SCAFF10X_SCAFF10X_MINSHAREDBARCODES} -ne 0 ]]
	then
		SCAFF10X_SCAFF10X_OPT="${SCAFF10X_SCAFF10X_OPT} -link ${SCAF_SCAFF10X_SCAFF10X_MINSHAREDBARCODES}"		
	fi
	
	if [[ -n ${SCAF_SCAFF10X_SCAFF10X_BLOCK} && ${SCAF_SCAFF10X_SCAFF10X_BLOCK} -ne 0 ]]
	then
		SCAFF10X_SCAFF10X_OPT="${SCAFF10X_SCAFF10X_OPT} -block ${SCAF_SCAFF10X_SCAFF10X_BLOCK}"		
	fi
	
    ### check input variable variables, and overwrite default pipeline if required
    if [[ -n ${SCAF_SCAFF10X_SCAFF10X_SAM} && -f ${SCAF_SCAFF10X_SCAFF10X_SAM} ]]
    then
    	SCAFF10X_SCAFF10X_OPT="${SCAFF10X_SCAFF10X_OPT} -sam ${SCAF_SCAFF10X_SCAFF10X_SAM}"
    elif [[ -n ${SCAF_SCAFF10X_SCAFF10X_BAM} && -f ${SCAF_SCAFF10X_SCAFF10X_BAM} ]]
    then
    	SCAFF10X_SCAFF10X_OPT="${SCAFF10X_SCAFF10X_OPT} -bam ${SCAF_SCAFF10X_SCAFF10X_BAM}"	
	fi
}

myTypes=("01_scaff10Xprepare, 02_scaff10Xscaff10x")
if [[ ${SC_SCAFF10X_TYPE} -eq 0 ]]
then 
    ### 01_scaff10Xprepare
    if [[ ${currentStep} -eq 1 ]]
    then
        ### clean up plans 
        for x in $(ls scaff10x_01_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        if [[ ! -d "${SC_SCAFF10X_READS}" ]]
        then
        	(>&2 echo "ERROR - set SC_SCAFF10X_READS to proper 10x read directory")
        	exit 1
   		fi
   		
		numR1Files=0
		for x in ${SC_SCAFF10X_READS}/${PROJECT_ID}_S*_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz
		do
			if [[ -f ${x} ]]
			then	
				numR1Files=$((${numR1Files}+1))	
			fi
		done
		
		if [[ ${numR1Files} -eq 0 ]]
        then
        	(>&2 echo "ERROR - cannot read 10x R1 files with following pattern: ${SC_SCAFF10X_READS}/${PROJECT_ID}_S*_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz")
        	exit 1
   		fi
   		
   		numR2Files=0
		for x in ${SC_SCAFF10X_READS}/${PROJECT_ID}_S*_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz
		do
			if [[ -f ${x} ]]
			then	
				numR2Files=$((${numR2Files}+1))	
			fi
		done
		
		if [[ ${numR2Files} -eq 0 ]]
        then
        	(>&2 echo "ERROR - cannot read 10x R2 files with following pattern: ${SC_SCAFF10X_READS}/${PROJECT_ID}_S*_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz")
        	exit 1
   		fi
   		
   		if [[ ${numR1Files} -ne ${numR2Files} ]]
        then
        	(>&2 echo "ERROR - 10x R1 files ${numR1Files} does not match R2 files ${numR2Files}")
        	exit 1
   		fi
   		
   		if [[ ! -f ${SC_SCAFF10X_REF} ]]
   		then
   			(>&2 echo "ERROR - set SC_SCAFF10X_REF to proper reference fasta file")
        	exit 1	
   		fi
   		
   		echo "if [[ -d ${SC_SCAFF10X_OUTDIR}/scaff10x_${SC_SCAFF10X_RUNID} ]]; then mv ${SC_SCAFF10X_OUTDIR}/scaff10x_${SC_SCAFF10X_RUNID} ${SC_SCAFF10X_OUTDIR}/scaff10x_${SC_SCAFF10X_RUNID}_$(date '+%Y-%m-%d_%H-%M-%S'); fi && mkdir -p ${SC_SCAFF10X_OUTDIR}/scaff10x_${SC_SCAFF10X_RUNID}" > scaff10x_01_scaff10Xprepare_single_${CONT_DB}.${slurmID}.plan
   		echo "mkdir -p ${SC_SCAFF10X_OUTDIR}/scaff10x_${SC_SCAFF10X_RUNID}/ref" >> scaff10x_01_scaff10Xprepare_single_${CONT_DB}.${slurmID}.plan
   		echo "mkdir -p ${SC_SCAFF10X_OUTDIR}/scaff10x_${SC_SCAFF10X_RUNID}/reads" >> scaff10x_01_scaff10Xprepare_single_${CONT_DB}.${slurmID}.plan
   		echo "ln -s -r ${SC_SCAFF10X_REF} ${SC_SCAFF10X_OUTDIR}/scaff10x_${SC_SCAFF10X_RUNID}/ref" >> scaff10x_01_scaff10Xprepare_single_${CONT_DB}.${slurmID}.plan
   		
   		for r1 in ${SC_SCAFF10X_READS}/${PROJECT_ID}_S[0-9]_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz
		do
			id=$(dirname ${r1})
			f1=$(basename ${r1})
			f2=$(echo "${f1}" | sed -e "s:_R1_:_R2_:")
			
			echo "ln -s -f ${id}/${f1} ${SC_SCAFF10X_OUTDIR}/scaff10x_${SC_SCAFF10X_RUNID}/reads"
			echo "ln -s -f ${id}/${f2} ${SC_SCAFF10X_OUTDIR}/scaff10x_${SC_SCAFF10X_RUNID}/reads"										
			echo "echo \"q1=${SC_SCAFF10X_OUTDIR}/scaff10x_${SC_SCAFF10X_RUNID}/reads/${f1}\" >> ${SC_SCAFF10X_OUTDIR}/scaff10x_${SC_SCAFF10X_RUNID}/scaff10x_inputReads.txt"
			echo "echo \"q2=${SC_SCAFF10X_OUTDIR}/scaff10x_${SC_SCAFF10X_RUNID}/reads/${f2}\" >> ${SC_SCAFF10X_OUTDIR}/scaff10x_${SC_SCAFF10X_RUNID}/scaff10x_inputReads.txt"							 
		done >> scaff10x_01_scaff10Xprepare_single_${CONT_DB}.${slurmID}.plan
		
		options="-debug 1 -tmp $(pwd)/${SC_SCAFF10X_OUTDIR}/scaff10x_${SC_SCAFF10X_RUNID}/"
		echo "${SCAFF10X_PATH}/scaff_reads ${options} ${SC_SCAFF10X_OUTDIR}/scaff10x_${SC_SCAFF10X_RUNID}/scaff10x_inputReads.txt scaff10x_BC_1.fastq scaff10x_BC_2.fastq" >> scaff10x_01_scaff10Xprepare_single_${CONT_DB}.${slurmID}.plan
		echo "scaff_reads $(cat ${SCAFF10X_PATH}/version.txt)" > scaff10x_01_scaff10Xprepare_single_${CONT_DB}.${slurmID}.version
    ### 02_scaff10Xscaff10x		
	elif [[ ${currentStep} -eq 2 ]]
    then
    	### clean up plans 
        for x in $(ls scaff10x_02_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        setScaff10xOptions
                  	
    	# add reference 
    	infiles="ref/$(basename ${SC_SCAFF10X_REF})"
    	if [[ -n ${SCAF_SCAFF10X_SCAFF10X_READSBC1} && -f ${SCAF_SCAFF10X_SCAFF10X_READSBC1} && -n ${SCAF_SCAFF10X_SCAFF10X_READSBC2} && -f ${SCAF_SCAFF10X_SCAFF10X_READSBC2} ]]
    	then
    		infiles="${infiles} ${SCAF_SCAFF10X_SCAFF10X_READSBC1} ${SCAF_SCAFF10X_SCAFF10X_READSBC2}"
    	else
    		infiles="${infiles} scaff10x_BC_1.fastq scaff10x_BC_2.fastq"
    	fi
                
        options="-debug 1 -tmp $(pwd)/${SC_SCAFF10X_OUTDIR}/scaff10x_${SC_SCAFF10X_RUNID}/"
        echo "${SCAFF10X_PATH}/scaff10x${SCAFF10X_SCAFF10X_OPT} ${options} ${infiles} ${PROJECT_ID}_m${FIX_FILT_OUTDIR}_x.p.fasta" > scaff10x_02_scaff10Xscaff10x_single_${CONT_DB}.${slurmID}.plan
		echo "scaff10X $(cat ${SCAFF10X_PATH}/version.txt)" > scaff10x_02_scaff10Xscaff10x_single_${CONT_DB}.${slurmID}.version
    else
        (>&2 echo "step ${currentStep} in SC_SCAFF10X_TYPE ${SC_SCAFF10X_TYPE} not supported")
        (>&2 echo "valid steps are: ${myTypes[${SC_SCAFF10X_TYPE}]}")
        exit 1            
    fi    		
else
    (>&2 echo "unknown SC_SCAFF10X_TYPE ${SC_SCAFF10X_TYPE}")
    (>&2 echo "supported types")
    x=0; while [ $x -lt ${#myTypes[*]} ]; do (>&2 echo "${myTypes[${x}]}"); done 
    exit 1
fi

exit 0