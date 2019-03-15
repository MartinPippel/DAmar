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

function setFastpOptions()
{
	if [[ -z "${CT_WHATSHAP_FASTP_THREADS}" ]]
	then 
		CT_WHATSHAP_FASTP_THREADS=4	
	fi
}

function setbwaOptions()
{
	CONTIG_BWA_OPT=""
	
	if [[ -z "${CT_WHATSHAP_BWA_THREADS}" ]]
	then 
		CT_WHATSHAP_BWA_THREADS=4	
	fi
	
	CONTIG_BWA_OPT="${CONTIG_BWA_OPT} -t ${CT_WHATSHAP_BWA_THREADS}"
	
	if [[ -n ${CT_WHATSHAP_BWA_VERBOSITY} ]]
	then 
		CONTIG_BWA_OPT="${CONTIG_BWA_OPT} -v ${CT_WHATSHAP_BWA_VERBOSITY}"
	fi
}


function setSamtoolsOptions()
{
	CONTIG_SAMTOOLS_OPT=""
	
	if [[ -z "${CT_WHATSHAP_SAMTOOLS_THREADS}" ]]
	then 
		CT_WHATSHAP_SAMTOOLS_THREADS=4	
	fi
	
	CONTIG_SAMTOOLS_OPT="${CONTIG_SAMTOOLS_OPT} -@ ${CT_WHATSHAP_SAMTOOLS_THREADS}"
	
	if [[ -n ${CT_WHATSHAP_SAMTOOLS_MEM} ]]
	then 
		CONTIG_SAMTOOLS_OPT="${CONTIG_SAMTOOLS_OPT} -m ${CT_WHATSHAP_SAMTOOLS_MEM}G"
	else
		CT_WHATSHAP_SAMTOOLS_MEM=2
		CONTIG_SAMTOOLS_OPT="${CONTIG_SAMTOOLS_OPT} -m ${CT_WHATSHAP_SAMTOOLS_MEM}G"
	fi
}

## type-1: Pacbio and 10x Genomcis reads 
myTypes=("1-WHprepareInput, 2-WHminimap2PacBio, 3-WHPacBioBamSplitByRef, 4-WHPacBioBamSeparate, 5-WHPacBioBamMerge")
if [[ ${CT_WHATSHAP_TYPE} -eq 0 ]]
then 
    ### 1-WHprepareInput
    if [[ ${currentStep} -eq 1 ]]
    then
        ### clean up plans 
        for x in $(ls whatshap_01_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        if [[ ! -f "${CT_WHATSHAP_REFFASTA}" ]]
        then
        	(>&2 echo "ERROR - set CT_WHATSHAP_REFFASTA to reference fasta file")
        	exit 1
   		fi
   		
   		echo "if [[ -d ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID} ]]; then mv ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID} ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}_$(date '+%Y-%m-%d_%H-%M-%S'); fi && mkdir ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}" > whatshap_01_WHprepareInput_block_${CONT_DB}.${slurmID}.plan
		echo "mkdir -p ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/pb_reads" >> whatshap_01_WHprepareInput_block_${CONT_DB}.${slurmID}.plan
		echo "mkdir -p ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/10x_reads" >> whatshap_01_WHprepareInput_block_${CONT_DB}.${slurmID}.plan
		echo "mkdir -p ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/pb_bams" >> whatshap_01_WHprepareInput_block_${CONT_DB}.${slurmID}.plan
		echo "mkdir -p ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/10x_bams" >> whatshap_01_WHprepareInput_block_${CONT_DB}.${slurmID}.plan
		echo "mkdir -p ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/ref && ln -s -r ${CT_WHATSHAP_REFFASTA} ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/ref && samtools faidx ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/ref/$(basename ${CT_WHATSHAP_REFFASTA}) && bwa index ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/ref/$(basename ${CT_WHATSHAP_REFFASTA})" >> whatshap_01_WHprepareInput_block_${CONT_DB}.${slurmID}.plan		
		
		# sanity checks
   		numFiles=0 
   		for x in ${CT_WHATSHAP_READS_PACBIO}/*.subreads.fa.gz   		
   		do
   			if [[ ! -f ${x} || ! -s ${x} ]]
   			then
   				(>&2 echo "WARNING - file ${x} not available or empty.")
   			else
   				numFiles=$((${numFiles}+1))
   				echo "ln -s -f -r ${x} ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/pb_reads/" >> whatshap_01_WHprepareInput_block_${CONT_DB}.${slurmID}.plan
   			fi      						
   		done

		echo "samtools $(${PACBIO_BASE_ENV} && samtools 2>&1 | grep Version | awk '{print $2}' && ${PACBIO_BASE_ENV_DEACT})" > whatshap_01_WHprepareInput_block_${CONT_DB}.${slurmID}.version
		echo "bwa $(${PACBIO_BASE_ENV} && bwa 2>&1 | grep Version | awk '{print $2}' && ${PACBIO_BASE_ENV_DEACT})" >> whatshap_01_WHprepareInput_block_${CONT_DB}.${slurmID}.version

		if [[ ${numFiles} -eq 0 ]]
		then
   		 
	   		for x in ${CT_WHATSHAP_READS_PACBIO}/*.subreads.bam   		
	   		do
	   			if [[ ! -f ${x} || ! -s ${x} ]]
	   			then
	   				(>&2 echo "WARNING - file ${x} not available or empty.")
	   			else
	   				numFiles=$((${numFiles}+1))
	   				# create dextract plans 
	   				echo "${DAZZLER_PATH}/bin/dextract -v -f -o ${x} | gzip > ${x%.subreads.bam}.fa.gz && ln -s -f -r ${x} ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/pb_reads/" >> whatshap_01_WHprepareInput_block_${CONT_DB}.${slurmID}.plan
	   				echo "DAZZLER $(git --git-dir=${DAZZLER_SOURCE_PATH}/DEXTRACTOR/.git rev-parse --short HEAD)" >> whatshap_01_WHprepareInput_block_${CONT_DB}.${slurmID}.version  	   				
	   			fi      						
	   		done
		fi
		
		if [[ ${numFiles} -eq 0 ]]
        then
        	(>&2 echo "ERROR - cannot find any pacbio sequel data with following pattern: ${CT_WHATSHAP_PACBIOFASTA}/*.subreads.bam")
        	exit 1
   		fi		
	### 2-WHPacBioMinimap2 
    elif [[ ${currentStep} -eq 2 ]]
    then
        ### clean up plans 
        for x in $(ls whatshap_02_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done

		if [[ ! -f ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/ref/$(basename ${CT_WHATSHAP_REFFASTA%.fasta}) ]]
        then
        	(>&2 echo "ERROR - Variable reference not found: ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/ref/$(basename ${CT_WHATSHAP_REFFASTA%.fasta})")
        	exit 1
        fi
       
        if [[ ! -f ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/ref/$(basename ${CT_WHATSHAP_REFFASTA}).fai ]]
        then
        	(>&2 echo "ERROR - Variable reference index not found: ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/ref/$(basename ${CT_WHATSHAP_REFFASTA}).fai")
        	exit 1
        fi
                        
        # sanity checks
   		numFiles=0 
   		for x in ${CT_WHATSHAP_PACBIOFASTA}/*.subreads.fa.gz   		
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
			(>&2 echo "ERROR - no input subreads.fa.gz file found")
	       	exit 1	
		fi
		
		ref="${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/ref/$(basename ${CT_WHATSHAP_REFFASTA%.fasta})"
				
        for x in ${CT_WHATSHAP_READS_PACBIO}/*.subreads.fa.gz   		
   		do
        	name=$(basename ${x%.subreads.fa.gz})
        	echo "minimap2 -a -x map-pb -R \"@RG\tID:${name}\tSM:${PROJECT_ID}_HIC\tLB:${PROJECT_ID}_HIC\tPL:PACBIO\tPU:none\" -t ${CT_WHATSHAP_MINIMAP2ALNTHREADS} ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/ref/${ref}.idx ${x} | samtools view -h - | samtools sort -@ ${CT_WHATSHAP_SAMTOOLSTHREADS} -m ${CT_WHATSHAP_SAMTOOLSMEM}G -o ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/pb_bams/${ref}_${name}_minimap2.sort.bam -T /tmp/${ref}_${name}_minimap2.sort.tmp && samtools index -@ ${CT_WHATSHAP_SAMTOOLSTHREADS} ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/pb_bams/${ref}_${name}_minimap2.sort.bam"        	
		done > whatshap_02_WHPacBioMinimap2_block_${CONT_DB}.${slurmID}.plan 
    	echo "minimap2 $(${PACBIO_BASE_ENV} && minimap2 --version && ${PACBIO_BASE_ENV_DEACT})" > whatshap_02_WHPacBioMinimap2_block_${CONT_DB}.${slurmID}.version
		echo "samtools $(${PACBIO_BASE_ENV} && samtools 2>&1 | grep Version | awk '{print $2}' && ${PACBIO_BASE_ENV_DEACT})" >> whatshap_02_WHPacBioMinimap2_block_${CONT_DB}.${slurmID}.version
	### 3-WHPacBioBamSplitByRef
    elif [[ ${currentStep} -eq 3 ]]
    then
    	### clean up plans 
        for x in $(ls whatshap_03_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        ### some sanity checks
        if [[ ! -f ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/ref/$(basename ${CT_WHATSHAP_REFFASTA}).fai ]]
        then
        	(>&2 echo "ERROR - Variable reference index not found: ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/ref/$(basename ${CT_WHATSHAP_REFFASTA%.fasta}).fai")
        	exit 1
        fi
        
        ref="${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/ref/$(basename ${CT_WHATSHAP_REFFASTA%.fasta})"
        
        for x in ${CT_WHATSHAP_READS_PACBIO}/*.subreads.fa.gz   		
   		do
        	name=$(basename ${x%.subreads.fa.gz})
        	if [[ ! -f ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/pb_bams/${ref}_${name}_minimap2.sort.bam ]]
        	then 
        		(>&2 echo "ERROR - missing file ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/pb_bams/${ref}_${name}_minimap2.sort.bam!")
        		(>&2 echo "        Does previous job WHPacBioMinimap2 finished properly?")
	       		exit 1	
        	fi        	
		done 
                		
		for x in ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/pb_bams/${ref}_*_minimap2.sort.bam   		
   		do
        	echo "bamtools split -in  ${x} -reference"			   			   		
		done > whatshap_03_WHPacBioBamSplitByRef_block_${CONT_DB}.${slurmID}.plan
		echo "$(${PACBIO_BASE_ENV} && bamtools --version | grep bamtools && ${PACBIO_BASE_ENV_DEACT})" > whatshap_03_WHPacBioBamSplitByRef_block_${CONT_DB}.${slurmID}.version
	### 4-WHPacBioBamSeparate 
    elif [[ ${currentStep} -eq 4 ]]
    then
        ### clean up plans 
        for x in $(ls whatshap_04_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        refIdx="${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/ref/$(basename ${CT_WHATSHAP_REFFASTA}).fai"
        
	    if [[ ! -f ${refIdx} ]]
	    then 
	    	(>&2 echo "ERROR - missing file ${refIdx}")
	       	exit 1
	    fi
	    
		echo "for x in \$(awk '{print \$1}' ${refIdx}); do mkdir -p ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/pb_bams/\$x; mv ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/pb_bams/*.REF_\${x}.bam ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/pb_bams/\${x}; find \$(pwd)/${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/pb_bams/\${x} -name \"*.bam\" > ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/pb_bams/\${x}/in.fof; done" > whatshap_04_WHPacBioBamSeparate_block_${CONT_DB}.${slurmID}.plan	    
	### 5-WHPacBioBamMerge 
    elif [[ ${currentStep} -eq 5 ]]
    then
        ### clean up plans 
        for x in $(ls whatshap_05_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        refIdx="${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/ref/$(basename ${CT_WHATSHAP_REFFASTA}).fai"
        
	    if [[ ! -f ${refIdx} ]]
	    then 
	    	(>&2 echo "ERROR - missing file ${refIdx}")
	       	exit 1
	    fi
	            	    
		# sanity checks
   		numFiles=0 
   		for x in $(awk '{print $1}' ${refIdx})   		
   		do
   			if [[ ! -f ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/pb_bams/${x}/in.fof || ! -s ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/pb_bams/${x}/in.fof ]]
   			then
   				(>&2 echo "WARNING - file ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/pb_bams/${x}/in.fof not available or empty")
   			else
   				numFiles=$((${numFiles}+1))
   			fi      						
   		done

		if [[ ${numFiles} -eq 0 ]]
		then
			(>&2 echo "ERROR - no splitted reference based alignment bam files found")
	       	exit 1	
		fi
		
		for x in $(awk '{print $1}' ${refIdx})   		
   		do
   			echo "bamtools merge -list ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/pb_bams/${x}/in.fof -out ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/pb_bams/${x}/ALL_${x}.bam"
		done > whatshap_05_WHPacBioBamMerge_block_${CONT_DB}.${slurmID}.plan
		echo "$(${PACBIO_BASE_ENV} && bamtools --version | grep bamtools && ${PACBIO_BASE_ENV_DEACT})" > whatshap_05_WHPacBioBamMerge_block_${CONT_DB}.${slurmID}.version    	
    else
        (>&2 echo "step ${currentStep} in CT_WHATSHAP_TYPE ${CT_WHATSHAP_TYPE} not supported")
        (>&2 echo "valid steps are: ${myTypes[${CT_WHATSHAP_TYPE}]}")
        exit 1            
    fi    		
else
    (>&2 echo "unknown CT_WHATSHAP_TYPE ${CT_WHATSHAP_TYPE}")
    (>&2 echo "supported types")
    x=0; while [ $x -lt ${#myTypes[*]} ]; do (>&2 echo "${myTypes[${x}]}"); done 
    exit 1
fi

exit 0