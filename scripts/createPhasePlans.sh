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
	if [[ -z "${CT_PHASE_FASTP_THREADS}" ]]
	then 
		CT_PHASE_FASTP_THREADS=4	
	fi
}

function setbwaOptions()
{
	CONTIG_BWA_OPT=""
	
	if [[ -z "${CT_PHASE_BWA_THREADS}" ]]
	then 
		CT_PHASE_BWA_THREADS=4	
	fi
	
	CONTIG_BWA_OPT="${CONTIG_BWA_OPT} -t ${CT_PHASE_BWA_THREADS}"
	
	if [[ -n ${CT_PHASE_BWA_VERBOSITY} ]]
	then 
		CONTIG_BWA_OPT="${CONTIG_BWA_OPT} -v ${CT_PHASE_BWA_VERBOSITY}"
	fi
}


function setSamtoolsOptions()
{
	CONTIG_SAMTOOLS_OPT=""
	
	if [[ -z "${CT_PHASE_SAMTOOLS_THREADS}" ]]
	then 
		CT_PHASE_SAMTOOLS_THREADS=4	
	fi
	
	CONTIG_SAMTOOLS_OPT="${CONTIG_SAMTOOLS_OPT} -@ ${CT_PHASE_SAMTOOLS_THREADS}"
	
	if [[ -n ${CT_PHASE_SAMTOOLS_MEM} ]]
	then 
		CONTIG_SAMTOOLS_OPT="${CONTIG_SAMTOOLS_OPT} -m ${CT_PHASE_SAMTOOLS_MEM}G"
	else
		CT_PHASE_SAMTOOLS_MEM=2
		CONTIG_SAMTOOLS_OPT="${CONTIG_SAMTOOLS_OPT} -m ${CT_PHASE_SAMTOOLS_MEM}G"
	fi
}

## type-1 [Whatshap]   - pacbio, 10x: 		01_WhatshapPrepareInput, 02_WhatshapMinimap2PacBio, 03_WhatshapPacBioBamSplitByRef, 04_WhatshapPacBioBamSplitByRef, 05_WhatshapPacBioBamMerge
## type-2 [Longranger] - 10x: 				01_LongrangerPrepareInput, 02_LongrangerLongrangerWgs, 03_LongrangerBcftoolsConsensus
## type-3 [HapCut2]    - pacbio, 10x, HiC: 	todo

myTypes=("01_WhatshapPrepareInput, 02_WhatshapMinimap2PacBio, 03_WhatshapPacBioBamSplitByRef, 04_WhatshapPacBioBamSplitByRef, 05_WhatshapPacBioBamMerge",
"01_LongrangerPrepareInput, 02_LongrangerLongrangerWgs, 03_LongrangerBcftoolsConsensus", "01_HapCut2PrepareInput")
if [[ ${CT_PHASE_TYPE} -eq 0 ]]
then 
    ### 01_WhatshapPrepareInput
    if [[ ${currentStep} -eq 1 ]]
    then
        ### clean up plans 
        for x in $(ls phase_01_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        if [[ ! -f "${CT_PHASE_REFFASTA}" ]]
        then
        	(>&2 echo "ERROR - set CT_PHASE_REFFASTA to reference fasta file")
        	exit 1
   		fi
   		
   		echo "if [[ -d ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID} ]]; then mv ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID} ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}_$(date '+%Y-%m-%d_%H-%M-%S'); fi && mkdir ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}" > phase_01_WhatshapPrepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "mkdir -p ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/pb_reads" >> phase_01_WhatshapPrepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "mkdir -p ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/10x_reads" >> phase_01_WhatshapPrepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "mkdir -p ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/pb_bams" >> phase_01_WhatshapPrepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "mkdir -p ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/10x_bams" >> phase_01_WhatshapPrepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "mkdir -p ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/ref && ln -s -r ${CT_PHASE_REFFASTA} ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/ref && samtools faidx ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/ref/$(basename ${CT_PHASE_REFFASTA}) && bwa index ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/ref/$(basename ${CT_PHASE_REFFASTA})" >> phase_01_WhatshapPrepareInput_single_${CONT_DB}.${slurmID}.plan		
		
		# sanity checks
   		numFiles=0 
   		for x in ${CT_PHASE_READS_PACBIO}/*.subreads.fa.gz   		
   		do
   			if [[ ! -f ${x} || ! -s ${x} ]]
   			then
   				(>&2 echo "WARNING - file ${x} not available or empty.")
   			else
   				numFiles=$((${numFiles}+1))
   				echo "ln -s -f -r ${x} ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/pb_reads/" >> phase_01_WhatshapPrepareInput_single_${CONT_DB}.${slurmID}.plan
   			fi      						
   		done

		echo "samtools $(${PACBIO_BASE_ENV} && samtools 2>&1 | grep Version | awk '{print $2}' && ${PACBIO_BASE_ENV_DEACT})" > phase_01_WhatshapPrepareInput_single_${CONT_DB}.${slurmID}.version
		echo "bwa $(${PACBIO_BASE_ENV} && bwa 2>&1 | grep Version | awk '{print $2}' && ${PACBIO_BASE_ENV_DEACT})" >> phase_01_WhatshapPrepareInput_single_${CONT_DB}.${slurmID}.version

		if [[ ${numFiles} -eq 0 ]]
		then
   		 
	   		for x in ${CT_PHASE_READS_PACBIO}/*.subreads.bam   		
	   		do
	   			if [[ ! -f ${x} || ! -s ${x} ]]
	   			then
	   				(>&2 echo "WARNING - file ${x} not available or empty.")
	   			else
	   				numFiles=$((${numFiles}+1))
	   				# create dextract plans 
	   				echo "${DAZZLER_PATH}/bin/dextract -v -f -o ${x} | gzip > ${x%.subreads.bam}.fa.gz && ln -s -f -r ${x} ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/pb_reads/" >> phase_01_WhatshapPrepareInput_single_${CONT_DB}.${slurmID}.plan
	   				echo "DAZZLER $(git --git-dir=${DAZZLER_SOURCE_PATH}/DEXTRACTOR/.git rev-parse --short HEAD)" >> phase_01_WhatshapPrepareInput_single_${CONT_DB}.${slurmID}.version  	   				
	   			fi      						
	   		done
		fi
		
		if [[ ${numFiles} -eq 0 ]]
        then
        	(>&2 echo "ERROR - cannot find any pacbio sequel data with following pattern: ${CT_PHASE_PACBIOFASTA}/*.subreads.bam")
        	exit 1
   		fi		
	### 02_WhatshapMinimap2PacBio 
    elif [[ ${currentStep} -eq 2 ]]
    then
        ### clean up plans 
        for x in $(ls phase_02_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done

		if [[ ! -f ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/ref/$(basename ${CT_PHASE_REFFASTA%.fasta}) ]]
        then
        	(>&2 echo "ERROR - Variable reference not found: ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/ref/$(basename ${CT_PHASE_REFFASTA%.fasta})")
        	exit 1
        fi
       
        if [[ ! -f ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/ref/$(basename ${CT_PHASE_REFFASTA}).fai ]]
        then
        	(>&2 echo "ERROR - Variable reference index not found: ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/ref/$(basename ${CT_PHASE_REFFASTA}).fai")
        	exit 1
        fi
                        
        # sanity checks
   		numFiles=0 
   		for x in ${CT_PHASE_PACBIOFASTA}/*.subreads.fa.gz   		
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
		
		ref="${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/ref/$(basename ${CT_PHASE_REFFASTA%.fasta})"
				
        for x in ${CT_PHASE_READS_PACBIO}/*.subreads.fa.gz   		
   		do
        	name=$(basename ${x%.subreads.fa.gz})
        	echo "minimap2 -a -x map-pb -R \"@RG\tID:${name}\tSM:${PROJECT_ID}_HIC\tLB:${PROJECT_ID}_HIC\tPL:PACBIO\tPU:none\" -t ${CT_PHASE_MINIMAP2ALNTHREADS} ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/ref/${ref}.idx ${x} | samtools view -h - | samtools sort -@ ${CT_PHASE_SAMTOOLSTHREADS} -m ${CT_PHASE_SAMTOOLSMEM}G -o ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/pb_bams/${ref}_${name}_minimap2.sort.bam -T /tmp/${ref}_${name}_minimap2.sort.tmp && samtools index -@ ${CT_PHASE_SAMTOOLSTHREADS} ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/pb_bams/${ref}_${name}_minimap2.sort.bam"        	
		done > phase_02_PhasePacBioMinimap2_block_${CONT_DB}.${slurmID}.plan 
    	echo "minimap2 $(${PACBIO_BASE_ENV} && minimap2 --version && ${PACBIO_BASE_ENV_DEACT})" > phase_02_PhasePacBioMinimap2_block_${CONT_DB}.${slurmID}.version
		echo "samtools $(${PACBIO_BASE_ENV} && samtools 2>&1 | grep Version | awk '{print $2}' && ${PACBIO_BASE_ENV_DEACT})" >> phase_02_PhasePacBioMinimap2_block_${CONT_DB}.${slurmID}.version
	### 03_WhatshapPacBioBamSplitByRef
    elif [[ ${currentStep} -eq 3 ]]
    then
    	### clean up plans 
        for x in $(ls phase_03_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        ### some sanity checks
        if [[ ! -f ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/ref/$(basename ${CT_PHASE_REFFASTA}).fai ]]
        then
        	(>&2 echo "ERROR - Variable reference index not found: ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/ref/$(basename ${CT_PHASE_REFFASTA%.fasta}).fai")
        	exit 1
        fi
        
        ref="${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/ref/$(basename ${CT_PHASE_REFFASTA%.fasta})"
        
        for x in ${CT_PHASE_READS_PACBIO}/*.subreads.fa.gz   		
   		do
        	name=$(basename ${x%.subreads.fa.gz})
        	if [[ ! -f ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/pb_bams/${ref}_${name}_minimap2.sort.bam ]]
        	then 
        		(>&2 echo "ERROR - missing file ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/pb_bams/${ref}_${name}_minimap2.sort.bam!")
        		(>&2 echo "        Does previous job 02_WhatshapMinimap2PacBio finished properly?")
	       		exit 1	
        	fi        	
		done 
                		
		for x in ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/pb_bams/${ref}_*_minimap2.sort.bam   		
   		do
        	echo "bamtools split -in  ${x} -reference"			   			   		
		done > phase_03_WhatshapPacBioBamSplitByRef_block_${CONT_DB}.${slurmID}.plan
		echo "$(${PACBIO_BASE_ENV} && bamtools --version | grep bamtools && ${PACBIO_BASE_ENV_DEACT})" > phase_03_WhatshapPacBioBamSplitByRef_block_${CONT_DB}.${slurmID}.version
	### 04_WhatshapPacBioBamSplitByRef 
    elif [[ ${currentStep} -eq 4 ]]
    then
        ### clean up plans 
        for x in $(ls phase_04_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        refIdx="${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/ref/$(basename ${CT_PHASE_REFFASTA}).fai"
        
	    if [[ ! -f ${refIdx} ]]
	    then 
	    	(>&2 echo "ERROR - missing file ${refIdx}")
	       	exit 1
	    fi
	    
		echo "for x in \$(awk '{print \$1}' ${refIdx}); do mkdir -p ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/pb_bams/\$x; mv ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/pb_bams/*.REF_\${x}.bam ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/pb_bams/\${x}; find \$(pwd)/${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/pb_bams/\${x} -name \"*.bam\" > ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/pb_bams/\${x}/in.fof; done" > phase_04_PhasePacBioBamSeparate_block_${CONT_DB}.${slurmID}.plan	    
	### 05_WhatshapPacBioBamMerge 
    elif [[ ${currentStep} -eq 5 ]]
    then
        ### clean up plans 
        for x in $(ls phase_05_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        refIdx="${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/ref/$(basename ${CT_PHASE_REFFASTA}).fai"
        
	    if [[ ! -f ${refIdx} ]]
	    then 
	    	(>&2 echo "ERROR - missing file ${refIdx}")
	       	exit 1
	    fi
	            	    
		# sanity checks
   		numFiles=0 
   		for x in $(awk '{print $1}' ${refIdx})   		
   		do
   			if [[ ! -f ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/pb_bams/${x}/in.fof || ! -s ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/pb_bams/${x}/in.fof ]]
   			then
   				(>&2 echo "WARNING - file ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/pb_bams/${x}/in.fof not available or empty")
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
   			echo "bamtools merge -list ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/pb_bams/${x}/in.fof -out ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/pb_bams/${x}/ALL_${x}.bam"
		done > phase_05_WhatshapPacBioBamMerge_block_${CONT_DB}.${slurmID}.plan
		echo "$(${PACBIO_BASE_ENV} && bamtools --version | grep bamtools && ${PACBIO_BASE_ENV_DEACT})" > phase_05_WhatshapPacBioBamMerge_block_${CONT_DB}.${slurmID}.version    	
    else
        (>&2 echo "step ${currentStep} in CT_PHASE_TYPE ${CT_PHASE_TYPE} not supported")
        (>&2 echo "valid steps are: ${myTypes[${CT_PHASE_TYPE}]}")
        exit 1            
    fi  
elif [[ ${CT_PHASE_TYPE} -eq 1 ]]
then 
    ### 01_LongrangerPrepareInput
    if [[ ${currentStep} -eq 1 ]]
    then
        ### clean up plans 
        for x in $(ls phase_01_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        if [[ ! -f "${CT_PHASE_REFFASTA}" ]]
        then
        	(>&2 echo "ERROR - set CT_PHASE_REFFASTA to reference fasta file")
        	exit 1
   		fi
   		
   		echo "if [[ -d ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID} ]]; then mv ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID} ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}_$(date '+%Y-%m-%d_%H-%M-%S'); fi && mkdir ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}" > phase_01_LongrangerPrepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "mkdir -p ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/bams" >> phase_01_LongrangerPrepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "mkdir -p ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/ref" >> phase_01_LongrangerPrepareInput_single_${CONT_DB}.${slurmID}.plan
		
		### prepare fasta file, i.e. sort according length and use only first 500 contigs for phasing
	
		NCONTIGS=500
		IN=${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/ref/phase.fasta
		IGNORE=${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/ref/ignore.fasta
		echo "awk '/^>/ {printf(\"%s%s\\t\",(N>0?\"\\n\":\"\"),\$0);N++;next;} {printf(\"%s\",\$0);} END {printf(\"\\n\");}' ${CT_PHASE_REFFASTA} | awk -F '\t' '{printf(\"%d\\t%s\\n\",length(\$2),\$0);}' | sort -k1,1rn | cut -f 2- | tee >(head -n ${NCONTIGS} | tr \"\\t\" \"\\n\" > ${IN}) | tail -n +${NCONTIGS} | tr \"\\t\" \"\\n\" > ${IGNORE}" >> phase_01_LongrangerPrepareInput_single_${CONT_DB}.${slurmID}.plan 
				
	### 02_LongrangerLongrangerWgs
    elif [[ ${currentStep} -eq 2 ]]
    then
        ### clean up plans 
        for x in $(ls phase_02_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        REFNAME=phase.fasta
        
        if [[ ! -f "${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/ref/${REFNAME}" ]]
        then
    		(>&2 echo "ERROR - cannot find reference fasta file \"${REFNAME}\" in dir \"${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/ref\"")
        	exit 1
   		fi
        
		if [[ ! -d ${TENX_PATH} ]]
        then 
        	(>&2 echo "ERROR - cannot find 10x reads. Variable TENX_PATH has to be set to a directoty containing 10x reads.")
        	exit 1
    	fi
    	 
        if [[ ! -d ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/ref/refdata-${REFNAME%.fasta} ]]
        then
        	echo "cd ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/ref && ${LONGRANGER_PATH}/longranger mkref ${REFNAME} && cd ../../../ " 
        	echo "cd ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/bams && ${LONGRANGER_PATH}/longranger wgs --id=10x_${PROJECT_ID}_longrangerWgs --fastqs=${TENX_PATH} --sample=${PROJECT_ID} --reference=../ref/refdata-${REFNAME%.fasta} --jobmode=slurm --localcores=24 --localmem=128 --maxjobs=500 --jobinterval=5000 --disable-ui --nopreflight && cd ../../../"
    	else 
    		(>&2 echo "[WARNING] Using previously created reference file ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/ref/refdata-${REFNAME}. Please remove that folder to rerun longranger mkref")
    		echo "cd ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/bams && ${LONGRANGER_PATH}/longranger wgs --vcmode=gatk:${GATK_PATH} --id=10x_${PROJECT_ID}_longrangerWgs --fastqs=${TENX_PATH} --sample=${PROJECT_ID} --reference=../ref/refdata-${REFNAME%.fasta} --jobmode=slurm --localcores=24 --localmem=128 --maxjobs=500 --jobinterval=5000 --disable-ui --nopreflight && cd ../../../"
    	fi > phase_02_LongrangerLongrangerWgs_single_${CONT_DB}.${slurmID}.plan                
        
        echo "$(${LONGRANGER_PATH}/longranger mkref --version | head -n1)" > phase_02_LongrangerLongrangerWgs_single_${CONT_DB}.${slurmID}.version
        echo "$(${LONGRANGER_PATH}/longranger wgs --version | head -n1)" >> phase_02_LongrangerLongrangerWgs_single_${CONT_DB}.${slurmID}.version   		
   		
	else
        (>&2 echo "step ${currentStep} in CT_PHASE_TYPE ${CT_PHASE_TYPE} not supported")
        (>&2 echo "valid steps are: ${myTypes[${CT_PHASE_TYPE}]}")
        exit 1            
    fi  
elif [[ ${CT_PHASE_TYPE} -eq 2 ]]
then 
    ### 01_HapCut2PrepareInput
    if [[ ${currentStep} -eq 1 ]]
    then
        ### clean up plans 
        for x in $(ls phase_01_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        if [[ ! -f "${CT_PHASE_REFFASTA}" ]]
        then
        	(>&2 echo "ERROR - set CT_PHASE_REFFASTA to reference fasta file")
        	exit 1
   		fi
	else
        (>&2 echo "step ${currentStep} in CT_PHASE_TYPE ${CT_PHASE_TYPE} not supported")
        (>&2 echo "valid steps are: ${myTypes[${CT_PHASE_TYPE}]}")
        exit 1            
    fi      
else
    (>&2 echo "unknown CT_PHASE_TYPE ${CT_PHASE_TYPE}")
    (>&2 echo "supported types")
    x=0; while [ $x -lt ${#myTypes[*]} ]; do (>&2 echo "${myTypes[${x}]}"); done 
    exit 1
fi

exit 0