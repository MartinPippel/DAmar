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

myTypes=("1-prepInFasta, 2-createMinimap2RefIndex, 3-minimap2, 4-bamMerge, 5-readCovHist, 6-contigCovHist, 7-purgeHaplotigs, 8-statistics")
if [[ ${CT_PURGEHAPLOTIGS_TYPE} -eq 0 ]]
then 
    ### 1-prepInFasta
    if [[ ${currentStep} -eq 1 ]]
    then
        ### clean up plans 
        for x in $(ls purgeHaplotigs_01_*_*_${FIX_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        if [[ ! -f "${CT_PURGEHAPLOTIGS_INFASTA}" ]]
        then
        	(>&2 echo "ERROR - set CT_PURGEHAPLOTIGS_INFASTA to input fasta file")
        	exit 1
   		fi
   		
   		echo "if [[ -d ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID} ]]; then mv ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID} ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}_$(date '+%Y-%m-%d_%H-%M-%S'); fi && mkdir ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}" > purgeHaplotigs_01_prepInFasta_single_${FIX_DB}.${slurmID}.plan
		echo "ln -s ${CT_PURGEHAPLOTIGS_INFASTA} ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}" >> purgeHaplotigs_01_prepInFasta_single_${FIX_DB}.${slurmID}.plan
	### 2-createMinimap2RefIndex
    elif [[ ${currentStep} -eq 2 ]]
    then
    	### clean up plans 
        for x in $(ls purgeHaplotigs_02_*_*_${FIX_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        ref=$(basename ${CT_PURGEHAPLOTIGS_INFASTA%.fasta})
        
        echo "minimap2 -t ${CT_PURGEHAPLOTIGS_MINIMAP2IDXTHREADS} -x map-pb -d ${ref}.idx ${ref}.fasta" > purgeHaplotigs_02_createMinimap2RefIndex_single_${FIX_DB}.${slurmID}.plan
    ### 3-minimap2
    elif [[ ${currentStep} -eq 3 ]]
    then
        ### clean up plans 
        for x in $(ls purgeHaplotigs_03_*_*_${FIX_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        if [[ ! -d ${CT_PURGEHAPLOTIGS_PACBIOFASTA} ]]
        then
        	(>&2 echo "ERROR - Variable ${CT_PURGEHAPLOTIGS_PACBIOFASTA} is not set or cannot be accessed")
        	exit 1
        fi
                
        # sanity checks
   		numFiles=0 
   		for x in ${CT_PURGEHAPLOTIGS_BAM}/*.subreads.fasta.gz   		
   		do
   			if [[ ! -f ${x} || ! -s ${x} ]]
   			then
   				(>&2 echo "WARNING - file ${x} not available or empty")
   			else
   				numFiles=$((${numFiles}+1))
   			fi      						
   		done

		if [[ ${numFiles} -eq 0 ]]
		then
			(>&2 echo "ERROR - no input bam file found")
	       	exit 1	
		fi
		
		ref=$(basename ${CT_PURGEHAPLOTIGS_INFASTA%.fasta})
		
		if [[ ! -f ${ref}.idx ]]
		then 
			(>&2 echo "ERROR - file ${ref}.idx not available. Create an reference index first!")
			exit 1
		fi
		
        for x in ${CT_PURGEHAPLOTIGS_BAM}/*.subreads.fasta.gz   		
   		do
        	name=$(basename ${x%.subreads.fasta.gz})
        	    		
    	echo "minimap2 -a -x map-pb -t ${CT_PURGEHAPLOTIGS_MINIMAP2ALNTHREADS} ${ref}.idx ${x} | samtools view -hF 256 - | samtools sort -@ ${CT_PURGEHAPLOTIGS_SAMTOOLSTHREADS} -m ${CT_PURGEHAPLOTIGS_SAMTOOLSMEM}G -o ${ref}_${name}_minmap2.sort.bam -T ${ref}_${name}_minmap2.sort.tmp"        	
		done > purgeHaplotigs_03_minimap2_block_${FIX_DB}.${slurmID}.plan 
    	echo "minimap2 $(minimap2 --version)" > purgeHaplotigs_03_minimap2_block_${FIX_DB}.${slurmID}.version
		echo "samtools $(samtools 2>&1 | grep Version | awk '{print $2}')" >> purgeHaplotigs_03_minimap2_block_${FIX_DB}.${slurmID}.version		
    else
        (>&2 echo "step ${currentStep} in CT_PURGEHAPLOTIGS_TYPE ${CT_PURGEHAPLOTIGS_TYPE} not supported")
        (>&2 echo "valid steps are: ${myTypes[${CT_PURGEHAPLOTIGS_TYPE}]}")
        exit 1            
    fi    		
else
    (>&2 echo "unknown CT_PURGEHAPLOTIGS_TYPE ${CT_PURGEHAPLOTIGS_TYPE}")
    (>&2 echo "supported types")
    x=0; while [ $x -lt ${#myTypes[*]} ]; do (>&2 echo "${myTypes[${x}]}"); done 
    exit 1
fi

exit 0