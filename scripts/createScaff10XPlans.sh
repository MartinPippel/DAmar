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

if [[ -z ${QUAST_PATH} || ! -f ${QUAST_PATH}/quast.py ]]
then
	(>&2 echo "Variable QUAST_PATH must be set to a proper quast installation directory!!")
    exit 1
fi

function setScaff10xOptions()
{
	SCAFF10X_SCAFF10X_OPT=""
	
	if [[ -n ${SC_SCAFF10X_SCAFF10X_THREADS} && ${SC_SCAFF10X_SCAFF10X_THREADS} -ne 0 ]]
	then
		SCAFF10X_SCAFF10X_OPT="${SCAFF10X_SCAFF10X_OPT} -nodes ${SC_SCAFF10X_SCAFF10X_THREADS}"		
	fi
	
	if [[ -n ${SC_SCAFF10X_SCAFF10X_ALIGNER} ]]
	then
		SCAFF10X_SCAFF10X_OPT="${SCAFF10X_SCAFF10X_OPT} -align ${SC_SCAFF10X_SCAFF10X_ALIGNER}"		
	fi
	
	if [[ -n ${SC_SCAFF10X_SCAFF10X_SCORE} && ${SC_SCAFF10X_SCAFF10X_SCORE} -ne 0 ]]
	then
		SCAFF10X_SCAFF10X_OPT="${SCAFF10X_SCAFF10X_OPT} -score ${SC_SCAFF10X_SCAFF10X_SCORE}"		
	fi
	
	if [[ -n ${SC_SCAFF10X_SCAFF10X_MATRIX} && ${SC_SCAFF10X_SCAFF10X_MATRIX} -ne 0 ]]
	then
		SCAFF10X_SCAFF10X_OPT="${SCAFF10X_SCAFF10X_OPT} -matrix ${SC_SCAFF10X_SCAFF10X_MATRIX}"		
	fi
	
	if [[ -n ${SC_SCAFF10X_SCAFF10X_MINREADS} && ${SC_SCAFF10X_SCAFF10X_MINREADS} -ne 0 ]]
	then
		SCAFF10X_SCAFF10X_OPT="${SCAFF10X_SCAFF10X_OPT} -reads ${SC_SCAFF10X_SCAFF10X_MINREADS}"		
	fi
	
	if [[ -n ${SC_SCAFF10X_SCAFF10X_LONGREAD} && ${SC_SCAFF10X_SCAFF10X_LONGREAD} -ne 0 ]]
	then
		SCAFF10X_SCAFF10X_OPT="${SCAFF10X_SCAFF10X_OPT} -longread ${SC_SCAFF10X_SCAFF10X_LONGREAD}"		
	fi
	
	if [[ -n ${SC_SCAFF10X_SCAFF10X_GAPSIZE} && ${SC_SCAFF10X_SCAFF10X_GAPSIZE} -ne 0 ]]
	then
		SCAFF10X_SCAFF10X_OPT="${SCAFF10X_SCAFF10X_OPT} -gap ${SC_SCAFF10X_SCAFF10X_GAPSIZE}"		
	fi
	
	if [[ -n ${SC_SCAFF10X_SCAFF10X_EDGELEN} && ${SC_SCAFF10X_SCAFF10X_EDGELEN} -ne 0 ]]
	then
		SCAFF10X_SCAFF10X_OPT="${SCAFF10X_SCAFF10X_OPT} -edge ${SC_SCAFF10X_SCAFF10X_EDGELEN}"		
	fi
	
	if [[ -n ${SC_SCAFF10X_SCAFF10X_MINSHAREDBARCODES} && ${SC_SCAFF10X_SCAFF10X_MINSHAREDBARCODES} -ne 0 ]]
	then
		SCAFF10X_SCAFF10X_OPT="${SCAFF10X_SCAFF10X_OPT} -link ${SC_SCAFF10X_SCAFF10X_MINSHAREDBARCODES}"		
	fi
	
	if [[ -n ${SC_SCAFF10X_SCAFF10X_BLOCK} && ${SC_SCAFF10X_SCAFF10X_BLOCK} -ne 0 ]]
	then
		SCAFF10X_SCAFF10X_OPT="${SCAFF10X_SCAFF10X_OPT} -block ${SC_SCAFF10X_SCAFF10X_BLOCK}"		
	fi
	
    ### check input variable variables, and overwrite default pipeline if required
    if [[ -n ${SC_SCAFF10X_SCAFF10X_SAM} && -f ${SC_SCAFF10X_SCAFF10X_SAM} ]]
    then
    	SCAFF10X_SCAFF10X_OPT="${SCAFF10X_SCAFF10X_OPT} -sam ${SC_SCAFF10X_SCAFF10X_SAM}"
    elif [[ -n ${SC_SCAFF10X_SCAFF10X_BAM} && -f ${SC_SCAFF10X_SCAFF10X_BAM} ]]
    then
    	SCAFF10X_SCAFF10X_OPT="${SCAFF10X_SCAFF10X_OPT} -bam ${SC_SCAFF10X_SCAFF10X_BAM}"	
	fi
}

function setBreak10xOptions()
{
	SCAFF10X_BREAK10X_OPT=""
	
	if [[ -n ${SC_SCAFF10X_BREAK10X_THREADS} && ${SC_SCAFF10X_BREAK10X_THREADS} -ne 0 ]]
	then
		SCAFF10X_BREAK10X_OPT="${SCAFF10X_BREAK10X_OPT} -nodes ${SC_SCAFF10X_BREAK10X_THREADS}"		
	fi
	
	if [[ -n ${SC_SCAFF10X_BREAK10X_READS} && ${SC_SCAFF10X_BREAK10X_READS} -ne 0 ]]
	then
		SCAFF10X_BREAK10X_OPT="${SCAFF10X_BREAK10X_OPT} -reads ${SC_SCAFF10X_BREAK10X_READS}"		
	fi

	if [[ -n ${SC_SCAFF10X_BREAK10X_SCORE} && ${SC_SCAFF10X_BREAK10X_SCORE} -ne 0 ]]
	then
		SCAFF10X_BREAK10X_OPT="${SCAFF10X_BREAK10X_OPT} -score ${SC_SCAFF10X_BREAK10X_SCORE}"		
	fi

	if [[ -n ${SC_SCAFF10X_BREAK10X_COVER} && ${SC_SCAFF10X_BREAK10X_COVER} -ne 0 ]]
	then
		SCAFF10X_BREAK10X_OPT="${SCAFF10X_BREAK10X_OPT} -cover ${SC_SCAFF10X_BREAK10X_COVER}"		
	fi
	
	if [[ -n ${SC_SCAFF10X_BREAK10X_RATIO} && ${SC_SCAFF10X_BREAK10X_RATIO} -ne 0 ]]
	then
		SCAFF10X_BREAK10X_OPT="${SCAFF10X_BREAK10X_OPT} -ratio ${SC_SCAFF10X_BREAK10X_RATIO}"		
	fi	
	
	if [[ -n ${SC_SCAFF10X_BREAK10X_GAP} && ${SC_SCAFF10X_BREAK10X_GAP} -ne 0 ]]
	then
		SCAFF10X_BREAK10X_OPT="${SCAFF10X_BREAK10X_OPT} -gap ${SC_SCAFF10X_BREAK10X_GAP}"		
	fi
}

myTypes=("01_scaff10Xprepare, 02_scaff10Xscaff10x, 03_scaff10Xbreak10x, 04_scaff10Xstatistics")
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
   		
   		if [[ -n ${SC_SCAFF10X_SCAFF10X_READSBC1} && -f ${SC_SCAFF10X_SCAFF10X_READSBC1} && -n ${SC_SCAFF10X_SCAFF10X_READSBC2} && -f ${SC_SCAFF10X_SCAFF10X_READSBC2} ]]
    	then
   			echo "using reads scaff10x_BC_1.fastq scaff10x_BC_2.fastq from a peevious run:"
   			echo "${SC_SCAFF10X_SCAFF10X_READSBC1}"
   			echo "${SC_SCAFF10X_SCAFF10X_READSBC2}"   			
   		else
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
		fi
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
    	# add BCR1 and BCR2 files
    	if [[ -n ${SC_SCAFF10X_SCAFF10X_READSBC1} && -f ${SC_SCAFF10X_SCAFF10X_READSBC1} && -n ${SC_SCAFF10X_SCAFF10X_READSBC2} && -f ${SC_SCAFF10X_SCAFF10X_READSBC2} ]]
    	then
    		### we need an absolute path if --tmp flag is used in scaff10x 
    		if [[ ! "${SC_SCAFF10X_SCAFF10X_READSBC1:0:1}" = "/" ]]
    		then 
    			SC_SCAFF10X_SCAFF10X_READSBC1=$(pwd)/${SC_SCAFF10X_SCAFF10X_READSBC1}
    		fi
    		
    		if [[ ! "${SC_SCAFF10X_SCAFF10X_READSBC2:0:1}" = "/" ]]
    		then 
    			SC_SCAFF10X_SCAFF10X_READSBC2=$(pwd)/${SC_SCAFF10X_SCAFF10X_READSBC2}
    		fi
    		infiles="${infiles} ${SC_SCAFF10X_SCAFF10X_READSBC1} ${SC_SCAFF10X_SCAFF10X_READSBC2}"
    	else
    		infiles="${infiles} scaff10x_BC_1.fastq scaff10x_BC_2.fastq"
    	fi
                
        options="-debug 1 -tmp $(pwd)/${SC_SCAFF10X_OUTDIR}/scaff10x_${SC_SCAFF10X_RUNID}/"
        echo "${SCAFF10X_PATH}/scaff10x${SCAFF10X_SCAFF10X_OPT} ${options} ${infiles} ${PROJECT_ID}_${SC_SCAFF10X_OUTDIR}_x.p.fasta" > scaff10x_02_scaff10Xscaff10x_single_${CONT_DB}.${slurmID}.plan
		echo "scaff10x $(cat ${SCAFF10X_PATH}/version.txt)" > scaff10x_02_scaff10Xscaff10x_single_${CONT_DB}.${slurmID}.version
	### 03_scaff10Xbreak10x		
	elif [[ ${currentStep} -eq 3 ]]
    then
		### clean up plans 
        for x in $(ls scaff10x_03_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        setBreak10xOptions
        # add reference 
    	inputAssembly="ref/$(basename ${SC_SCAFF10X_REF})"
    	inputScaffold="${PROJECT_ID}_${SC_SCAFF10X_OUTDIR}_x.p.fasta"
    	# add BCR1 and BCR2 files
    	if [[ -n ${SC_SCAFF10X_SCAFF10X_READSBC1} && -f ${SC_SCAFF10X_SCAFF10X_READSBC1} && -n ${SC_SCAFF10X_SCAFF10X_READSBC2} && -f ${SC_SCAFF10X_SCAFF10X_READSBC2} ]]
    	then
    		### we need an absolute path if --tmp flag is used in scaff10x 
    		if [[ ! "${SC_SCAFF10X_SCAFF10X_READSBC1:0:1}" = "/" ]]
    		then 
    			SC_SCAFF10X_SCAFF10X_READSBC1=$(pwd)/${SC_SCAFF10X_SCAFF10X_READSBC1}
    		fi
    		
    		if [[ ! "${SC_SCAFF10X_SCAFF10X_READSBC2:0:1}" = "/" ]]
    		then 
    			SC_SCAFF10X_SCAFF10X_READSBC2=$(pwd)/${SC_SCAFF10X_SCAFF10X_READSBC2}
    		fi
    		inputAssembly="${inputAssembly} ${SC_SCAFF10X_SCAFF10X_READSBC1} ${SC_SCAFF10X_SCAFF10X_READSBC2}"
    		inputScaffold="${inputScaffold} ${SC_SCAFF10X_SCAFF10X_READSBC1} ${SC_SCAFF10X_SCAFF10X_READSBC2}"
    	else
    		inputAssembly="${inputAssembly} scaff10x_BC_1.fastq scaff10x_BC_2.fastq"
    		inputScaffold="${inputScaffold} scaff10x_BC_1.fastq scaff10x_BC_2.fastq"
    	fi
                
        options="-debug 1 -tmp $(pwd)/${SC_SCAFF10X_OUTDIR}/scaff10x_${SC_SCAFF10X_RUNID}/"
        
		echo "${SCAFF10X_PATH}/scaff10x${SCAFF10X_BREAK10X_OPT} ${options} ${inputAssembly} $(basename ${SC_SCAFF10X_REF%.*})_break10x.fasta $(basename ${SC_SCAFF10X_REF%.*})_break10x.breaks" > scaff10x_03_scaff10Xbreak10x_block_${CONT_DB}.${slurmID}.plan
		echo "${SCAFF10X_PATH}/scaff10x${SCAFF10X_BREAK10X_OPT} ${options} ${inputScaffold} ${PROJECT_ID}_${SC_SCAFF10X_OUTDIR}_x.p_break10x.fasta ${PROJECT_ID}_${SC_SCAFF10X_OUTDIR}_x.p_break10x.breaks" >> scaff10x_03_scaff10Xbreak10x_block_${CONT_DB}.${slurmID}.plan
		echo "break10x $(cat ${SCAFF10X_PATH}/version.txt)" > scaff10x_03_scaff10Xbreak10x_block_${CONT_DB}.${slurmID}.version
	### 04_scaff10Xstatistics		
	elif [[ ${currentStep} -eq 4 ]]
    then
		### clean up plans 
        for x in $(ls scaff10x_04_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
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
    	echo "${SUBMIT_SCRIPTS_PATH}/assemblyStats.sh ${configFile} 12" > scaff10x_04_scaff10Xstatistics_single_${CONT_DB}.${slurmID}.plan
    	echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scaff10x_04_scaff10Xstatistics_single_${CONT_DB}.${slurmID}.version	
		echo "$(quast.py --version)" >> scaff10x_04_scaff10Xstatistics_single_${CONT_DB}.${slurmID}.version
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