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
	if [[ -z "${CT_FREEBAYES_FASTP_THREADS}" ]]
	then 
		CT_FREEBAYES_FASTP_THREADS=4	
	fi
}

function setbwaOptions()
{
	CONTIG_BWA_OPT=""
	
	if [[ -z "${CT_FREEBAYES_BWA_THREADS}" ]]
	then 
		CT_FREEBAYES_BWA_THREADS=4	
	fi
	
	CONTIG_BWA_OPT="${CONTIG_BWA_OPT} -t ${CT_FREEBAYES_BWA_THREADS}"
	
	if [[ -n ${CT_FREEBAYES_BWA_VERBOSITY} ]]
	then 
		CONTIG_BWA_OPT="${CONTIG_BWA_OPT} -v ${CT_FREEBAYES_BWA_VERBOSITY}"
	fi
}


function setSamtoolsOptions()
{
	CONTIG_SAMTOOLS_OPT=""
	
	if [[ -z "${CT_FREEBAYES_SAMTOOLS_THREADS}" ]]
	then 
		CT_FREEBAYES_SAMTOOLS_THREADS=4	
	fi
	
	CONTIG_SAMTOOLS_OPT="${CONTIG_SAMTOOLS_OPT} -@ ${CT_FREEBAYES_SAMTOOLS_THREADS}"
	
	if [[ -n ${CT_FREEBAYES_SAMTOOLS_MEM} ]]
	then 
		CONTIG_SAMTOOLS_OPT="${CONTIG_SAMTOOLS_OPT} -m ${CT_FREEBAYES_SAMTOOLS_MEM}G"
	else
		CT_FREEBAYES_SAMTOOLS_MEM=2
		CONTIG_SAMTOOLS_OPT="${CONTIG_SAMTOOLS_OPT} -m ${CT_FREEBAYES_SAMTOOLS_MEM}G"
	fi
}

function setPicardOptions()
{
	CONTIG_PICARD_OPT=""
	
	if [[ -z ${CT_FREEBAYES_PICARD_XMX} ]]
	then 
		CT_FREEBAYES_PICARD_XMX=8	
	fi
	
	if [[ -z ${CT_FREEBAYES_PICARD_XMS} ]]
	then 
		CT_FREEBAYES_PICARD_XMS=8	
	fi
	
	CONTIG_PICARD_OPT="-Xmx${CT_FREEBAYES_PICARD_XMX}G -Xms${CT_FREEBAYES_PICARD_XMS}G -XX:-UseGCOverheadLimit"	
}

if [[ -z ${LONGRANGER_PATH} || ! -f ${LONGRANGER_PATH}/longranger ]]
then
	(>&2 echo "Variable LONGRANGER_PATH must be set to a proper longranger installation directory!!")
    exit 1
fi

# Type: 0 [bwa mapping] - 01_FBprepareInput, 02_FBfastp, 03_FBbwa, 04_FBmarkDuplicates, 05_FBfreebayes, 06_FBconsensus, 07_FBstatistics 
# Type: 1 [longranger mapping] - 01_FBprepareInput, 02_FBlongrangerAlign, 03_FBfreebayes, 04_FBconsensus, 05_FBstatistics
myTypes=("01_FBprepareInput, 02_FBfastp, 03_FBbwa, 04_FBmarkDuplicates, 05_FBfreebayes, 06_FBconsensus, 07_FBstatistics", 
"01_FBprepareInput, 02_FBlongrangerAlign, 03_FBfreebayes, 04_FBconsensus, 05_FBstatistics")
if [[ ${CT_FREEBAYES_TYPE} -eq 0 ]]
then 
    ### 1-FBprepareInput
    if [[ ${currentStep} -eq 1 ]]
    then
        ### clean up plans 
        for x in $(ls freebayes_01_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        if [[ ! -f "${CT_FREEBAYES_REFFASTA}" ]]
        then
        	(>&2 echo "ERROR - set CT_FREEBAYES_REFFASTA to reference fasta file")
        	exit 1
   		fi
   		
   		echo "if [[ -d ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID} ]]; then mv ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID} ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}_$(date '+%Y-%m-%d_%H-%M-%S'); fi && mkdir ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}" > freebayes_01_FBprepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "mkdir -p ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/reads" >> freebayes_01_FBprepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "mkdir -p ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/bams" >> freebayes_01_FBprepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "mkdir -p ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/ref" >> freebayes_01_FBprepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "mkdir -p ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/freebayes" >> freebayes_01_FBprepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "ln -s -r ${CT_FREEBAYES_REFFASTA} ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/ref" >> freebayes_01_FBprepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "samtools faidx ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/ref/$(basename ${CT_FREEBAYES_REFFASTA})" >> freebayes_01_FBprepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "bwa index ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/ref/$(basename ${CT_FREEBAYES_REFFASTA})" >> freebayes_01_FBprepareInput_single_${CONT_DB}.${slurmID}.plan
		
		echo "samtools $(${CONDA_BASE_ENV} && samtools 2>&1 | grep Version | awk '{print $2}' && conda deactivate)" > freebayes_01_FBprepareInput_single_${CONT_DB}.${slurmID}.version
		echo "bwa $(${CONDA_BASE_ENV} && bwa 2>&1 | grep Version | awk '{print $2}' && conda deactivate)" >> freebayes_01_FBprepareInput_single_${CONT_DB}.${slurmID}.version
	### 2-FBfastp 
    elif [[ ${currentStep} -eq 2 ]]
    then
        ### clean up plans 
        for x in $(ls freebayes_02_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
		
		if [[ ! -d "${CT_FREEBAYES_READS}" ]]
        then
        	(>&2 echo "ERROR - set CT_FREEBAYES_READS to directory that contain the PROJECT_ID*.fastq.qz read files")
        	exit 1
   		fi
   		   				
		# we need to trim off barcode and adapter for 10x 
		if [[ "x${CT_FREEBAYES_READSTYPE}" == "x10x" ]]
		then
	   		numR1Files=0
			for x in ${CT_FREEBAYES_READS}/${PROJECT_ID}_S*_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz
			do
				if [[ -f ${x} ]]
				then	
					numR1Files=$((${numR1Files}+1))	
				fi
			done
			
			if [[ ${numR1Files} -eq 0 ]]
	        then
	        	(>&2 echo "ERROR - cannot read 10x R1 files with following pattern: ${CT_FREEBAYES_READS}/${PROJECT_ID}_S*_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz")
	        	exit 1
	   		fi
	   		
	   		numR2Files=0
			for x in ${CT_FREEBAYES_READS}/${PROJECT_ID}_S*_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz
			do
				if [[ -f ${x} ]]
				then	
					numR2Files=$((${numR2Files}+1))	
				fi
			done
			
			if [[ ${numR2Files} -eq 0 ]]
	        then
	        	(>&2 echo "ERROR - cannot read 10x R2 files with following pattern: ${CT_FREEBAYES_READS}/${PROJECT_ID}_S*_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz")
	        	exit 1
	   		fi
	   		
	   		if [[ ${numR1Files} -ne ${numR2Files} ]]
	        then
	        	(>&2 echo "ERROR - 10x R1 files ${numR1Files} does not match R2 files ${numR2Files}")
	        	exit 1
	   		fi
	   		
	   		setFastpOptions
	   		
			for r1 in ${CT_FREEBAYES_READS}/${PROJECT_ID}_S[0-9]_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz
			do
				id=$(dirname ${r1})
				f1=$(basename ${r1})
				f2=$(echo "${f1}" | sed -e "s:_R1_:_R2_:")
				o="${f1%_R1_???.fastq.gz}"											
				
				echo "fastp -i ${id}/${f1} -I ${id}/${f2} -f 23 -G -Q -j ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/reads/${o}.json -h ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/reads/${o}.html -w ${CT_FREEBAYES_FASTP_THREADS} -o ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/reads/${f1} -O ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/reads/${f2}"				 
			done > freebayes_02_FBfastp_block_${CONT_DB}.${slurmID}.plan
	   		echo "fastp $(${CONDA_BASE_ENV} && fastp 2>&1 | grep version | awk '{print $2}' && conda deactivate)" > freebayes_02_FBfastp_block_${CONT_DB}.${slurmID}.version	   		
		else		
			(>&2 echo "ERROR - set CT_FREEBAYES_READSTYPE ${CT_FREEBAYES_READSTYPE} not supported yet")
	        exit 1			
		fi
	### 3-FBbwa
    elif [[ ${currentStep} -eq 3 ]]
    then
    	### clean up plans 
        for x in $(ls freebayes_03_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        if [[ ! -d "${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/reads" ]]
        then
        	(>&2 echo "ERROR - cannot access directory ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/reads!")
        	exit 1
   		fi
   		
   		ref=${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/ref/$(basename ${CT_FREEBAYES_REFFASTA})
   		
   		if [[ ! -f "${ref}" ]]
        then
        (>&2 echo "ERROR - cannot reference fasta file ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/ref/$(basename ${CT_FREEBAYES_REFFASTA})!")
        	exit 1
   		fi
   		   				
		# we need to trim off barcode and adapter for 10x 
		if [[ "x${CT_FREEBAYES_READSTYPE}" == "x10x" ]]
		then
	   		numR1Files=0
			for x in ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/reads/${PROJECT_ID}_S*_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz
			do
				if [[ -f ${x} ]]
				then	
					numR1Files=$((${numR1Files}+1))	
				fi
			done
			
			if [[ ${numR1Files} -eq 0 ]]
	        then
	        	(>&2 echo "ERROR - cannot read 10x R1 files with following pattern: ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/reads/${PROJECT_ID}_S*_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz")
	        	exit 1
	   		fi
	   		
	   		numR2Files=0
			for x in ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/reads/${PROJECT_ID}_S*_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz
			do
				if [[ -f ${x} ]]
				then	
					numR2Files=$((${numR2Files}+1))	
				fi
			done
			
			if [[ ${numR2Files} -eq 0 ]]
	        then
	        	(>&2 echo "ERROR - cannot read 10x R2 files with following pattern: ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/reads/${PROJECT_ID}_S*_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz")
	        	exit 1
	   		fi
	   		
	   		if [[ ${numR1Files} -ne ${numR2Files} ]]
	        then
	        	(>&2 echo "ERROR - 10x R1 files ${numR1Files} does not match R2 files ${numR2Files}")
	        	exit 1
	   		fi
	   		
	   		setbwaOptions
	   		setSamtoolsOptions
	   		
			for r1 in ${CT_FREEBAYES_READS}/${PROJECT_ID}_S[0-9]_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz
			do
				id=$(dirname ${r1})
				f1=$(basename ${r1})
				f2=$(echo "${f1}" | sed -e "s:_R1_:_R2_:")
				o="${f1%_R1_???.fastq.gz}"											
				
				echo "bwa mem${CONTIG_BWA_OPT} -R \"@RG\tID:${PROJECT_ID}\tSM:10X\" ${ref} ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/reads/${f1} ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/reads/${f2} | samtools sort${CONTIG_SAMTOOLS_OPT} -O BAM -o ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/bams/${o}_bwa.bam && samtools index -@ ${CT_FREEBAYES_SAMTOOLS_THREADS} ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/bams/${o}_bwa.bam" 				 
			done > freebayes_03_FBbwa_block_${CONT_DB}.${slurmID}.plan
			
	   		echo "bwa $(${CONDA_BASE_ENV} && bwa 2>&1 | grep Version | awk '{print $2}' && conda deactivate)" > freebayes_03_FBbwa_block_${CONT_DB}.${slurmID}.version
        	echo "samtools $(${CONDA_BASE_ENV} && samtools 2>&1 | grep Version | awk '{print $2}' && conda deactivate)" >> freebayes_03_FBbwa_block_${CONT_DB}.${slurmID}.version
		else		
			(>&2 echo "ERROR - set CT_FREEBAYES_READSTYPE ${CT_FREEBAYES_READSTYPE} not supported yet")
	        exit 1			
		fi        	
    ### 4-FBmarkDuplicates
    elif [[ ${currentStep} -eq 4 ]]
    then
        ### clean up plans 
        for x in $(ls freebayes_04_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done               

		if [[ ! -d "${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/bams" ]]
        then
        	(>&2 echo "ERROR - cannot access directory ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/bams!")
        	exit 1
   		fi
   		
   		## if multiple bam files are available (e.g. different Lanes) then merge files prior to markduplicates
   		files=$(ls ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/bams/*_bwa.bam)
   		
   		setPicardOptions
   		
   		echo "picard MarkDuplicates $(${CONDA_BASE_ENV} && picard MarkDuplicates --version && conda deactivate)" > freebayes_04_FBmarkDuplicates_block_${CONT_DB}.${slurmID}.version
   		if [[ $(echo $files | wc -w) -eq 1 ]]
   		then
   			ob="${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/bams/${PROJECT_ID}_final10x.bam"
			m="${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/bams/${PROJECT_ID}_final10x.metrics"
   			
   			echo "picard ${CONTIG_PICARD_OPT} MarkDuplicates I=${files} O=${ob} M=${m} && samtools index -@ ${CT_FREEBAYES_SAMTOOLS_THREADS} ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/bams/${ob}" 				 
   		elif [[ $(echo $files | wc -w) -gt 1 ]]
   		then
   			mrg=${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/bams/${PROJECT_ID}_merged10x.bam
   			o=${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/bams/${PROJECT_ID}_final10x.bam
   			m=${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/bams/${PROJECT_ID}_final10x.metrics
   			i=$(echo -e ${files} | sed -e "s:${CT_FREEBAYES_OUTDIR}:I=${CT_FREEBAYES_OUTDIR}:g")
   			echo "picard ${CONTIG_PICARD_OPT} MergeSamFiles ${i} OUTPUT=${mrg} USE_THREADING=TRUE ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT && picard ${CONTIG_PICARD_OPT} MarkDuplicates I=${mrg} O=${o} M=${m} && samtools index -@ ${CT_FREEBAYES_SAMTOOLS_THREADS} ${o}"
   			echo "picard MergeSamFiles $(${CONDA_BASE_ENV} && picard MarkDuplicates --version && conda deactivate)" >> freebayes_04_FBmarkDuplicates_block_${CONT_DB}.${slurmID}.version	   			
		else
   	 		(>&2 echo "ERROR - cannot find file with following pattern: ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/bams/*_bwa.bam")
        	exit 1
   	 	fi > freebayes_04_FBmarkDuplicates_block_${CONT_DB}.${slurmID}.plan   			   		   		   						   	
	### 5-FBfreebayes
    elif [[ ${currentStep} -eq 5 ]]
    then
        ### clean up plans 
        for x in $(ls freebayes_05_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        if [[ ! -d "${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/bams" ]]
        then
        	(>&2 echo "ERROR - cannot access directory ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/bams!")
        	exit 1
   		fi
   		
   		outdir="${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/freebayes/"
   		if [[ ! -d "${outdir}" ]]
        then
        	(>&2 echo "ERROR - cannot access directory ${outdir}!")
        	exit 1
   		fi   		
        
        ref=${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/ref/$(basename ${CT_FREEBAYES_REFFASTA})
   		
   		if [[ ! -f "${ref}.fai" ]]
        then
        	(>&2 echo "ERROR - cannot reference fasta index file ${ref}.fai!")
        	exit 1
   		fi
   		
   		bam=${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/bams/${PROJECT_ID}_final10x.bam
   		if [[ ! -f ${bam} ]]
   		then 
   			(>&2 echo "ERROR - cannot final duplicate marked bam file: ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/bams/${PROJECT_ID}_final10x.bam!")
        	exit 1
   		fi 
   		
   		echo "$(awk -v bam=${bam} -v ref=${ref} -v out=${outdir} '{print "freebayes --bam "bam" --region "$1":1-"$2" -f "ref" | bcftools view --no-version -Ob -o "out$1":1-"$2".bcf"}' ${ref}.fai)" > freebayes_05_FBfreebayes_block_${CONT_DB}.${slurmID}.plan

		echo "freebayes $(${CONDA_BASE_ENV} && freebayes --version && conda deactivate)" > freebayes_05_FBfreebayes_block_${CONT_DB}.${slurmID}.version
		echo "bcftools $(${CONDA_BASE_ENV} && bcftools --version | head -n1 | awk '{print $2}' && conda deactivate)" >> freebayes_05_FBfreebayes_block_${CONT_DB}.${slurmID}.version
   	### 06_FBconsensus
    elif [[ ${currentStep} -eq 6 ]]
    then
        ### clean up plans 
        for x in $(ls freebayes_06_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
		ref=${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/ref/$(basename ${CT_FREEBAYES_REFFASTA})
   		
   		if [[ ! -f "${ref}.fai" ]]
        then
        	(>&2 echo "ERROR - cannot reference fasta index file ${ref}.fai!")
        	exit 1
   		fi
   		
   		indir="${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/freebayes/"
   		if [[ ! -d "${indir}" ]]
        then
        	(>&2 echo "ERROR - cannot access directory ${indir}!")
        	exit 1
   		fi
   		
   		outdir="${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/"
        
    	# create list of bcf files, same order as in ref.fai
        echo "awk -v d=\"${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/freebayes/\" '{print d\$1\":1-\"\$2\".bcf\"}' ${ref}.fai > ${outdir}/${PROJECT_ID}_10x_concatList.txt" > freebayes_06_FBconsensus_single_${CONT_DB}.${slurmID}.plan
        echo "bcftools concat -nf ${outdir}/${PROJECT_ID}_10x_concatList.txt | bcftools view -Ou -e'type=\"ref\"' | bcftools norm -Ob -f $ref -o ${outdir}/${PROJECT_ID}_10x.bcf" >> freebayes_06_FBconsensus_single_${CONT_DB}.${slurmID}.plan
        echo "bcftools index ${outdir}/${PROJECT_ID}_10x.bcf" >> freebayes_06_FBconsensus_single_${CONT_DB}.${slurmID}.plan
        echo "bcftools consensus -i'QUAL>1 && (GT=\"AA\" || GT=\"Aa\")' -Hla -f ${ref} ${outdir}/${PROJECT_ID}_10x.bcf > ${outdir}/${PROJECT_ID}_10x.fasta" >> freebayes_06_FBconsensus_single_${CONT_DB}.${slurmID}.plan
      
      	echo "bcftools $(${CONDA_BASE_ENV} && bcftools --version | head -n1 | awk '{print $2}' && conda deactivate)" > freebayes_06_FBconsensus_single_${CONT_DB}.${slurmID}.version  
    ### 7-FBstatistics
    elif [[ ${currentStep} -eq 7 ]]
    then
        ### clean up plans 
        for x in $(ls freebayes_07_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
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
    	echo "${SUBMIT_SCRIPTS_PATH}/assemblyStats.sh ${configFile} 11" > freebayes_07_FBstatistics_single_${CONT_DB}.${slurmID}.plan
    	git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD > freebayes_07_FBstatistics_single_${CONT_DB}.${slurmID}.version
    else
        (>&2 echo "step ${currentStep} in CT_FREEBAYES_TYPE ${CT_FREEBAYES_TYPE} not supported")
        (>&2 echo "valid steps are: ${myTypes[${CT_FREEBAYES_TYPE}]}")
        exit 1            
    fi    
elif [[ ${CT_FREEBAYES_TYPE} -eq 1 ]]
then 
    ### 01_FBprepareInput
    if [[ ${currentStep} -eq 1 ]]
    then
        ### clean up plans 
        for x in $(ls freebayes_01_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        if [[ ! -f "${CT_FREEBAYES_REFFASTA}" ]]
        then
        	(>&2 echo "ERROR - set CT_FREEBAYES_REFFASTA to reference fasta file")
        	exit 1
   		fi
   		
   		echo "if [[ -d ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID} ]]; then mv ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID} ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}_$(date '+%Y-%m-%d_%H-%M-%S'); fi && mkdir ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}" > freebayes_01_FBprepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "mkdir -p ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/bams" >> freebayes_01_FBprepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "mkdir -p ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/ref" >> freebayes_01_FBprepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "mkdir -p ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/freebayes" >> freebayes_01_FBprepareInput_single_${CONT_DB}.${slurmID}.plan
		# get rid of any colon's, as those will cause a crash of longranger		
		echo "sed -e \"s/:/-/g\" ${CT_FREEBAYES_REFFASTA} > ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/ref/$(basename ${CT_FREEBAYES_REFFASTA})" >> freebayes_01_FBprepareInput_single_${CONT_DB}.${slurmID}.plan		
	### 02_FBlongrangerAlign
    elif [[ ${currentStep} -eq 2 ]]
    then
        ### clean up plans 
        for x in $(ls freebayes_02_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        REFNAME=$(basename ${CT_FREEBAYES_REFFASTA})
        
        if [[ ! -f "${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/ref/${REFNAME}" ]]
        then
    		(>&2 echo "ERROR - cannot find reference fasta file \"${REFNAME}\" in dir \"${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/ref\"")
        	exit 1
   		fi
        
		if [[ ! -d ${TENX_PATH} ]]
        then 
        	(>&2 echo "ERROR - cannot find 10x reads. Variable TENX_PATH has to be set to a directoty containing 10x reads.")
        	exit 1
    	fi
    	 
        if [[ ! -d ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/ref/refdata-${REFNAME%.fasta} ]]
        then
        	echo "cd ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/ref && ${LONGRANGER_PATH}/longranger mkref ${REFNAME} && cd ../../../ " 
        	echo "cd ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/bams && ${LONGRANGER_PATH}/longranger align --id=10x_${PROJECT_ID}_longrangerAlign --fastqs=${TENX_PATH} --sample=${PROJECT_ID} --reference=../ref/refdata-${REFNAME%.fasta} --jobmode=slurm --localcores=38 --localmem=128 --maxjobs=1000 --jobinterval=5000 --disable-ui --nopreflight && cd ../../../"
    	else 
    		(>&2 echo "[WARNING] Using previously created reference file ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/ref/refdata-${REFNAME}. Please remove that folder to rerun longranger mkref" )
    		echo "cd ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/bams && ${LONGRANGER_PATH}/longranger align --id=10x_${PROJECT_ID}_longrangerAlign --fastqs=${TENX_PATH} --sample=${PROJECT_ID} --reference=../ref/refdata-${REFNAME%.fasta} --jobmode=slurm --localcores=38 --localmem=128 --maxjobs=1000 --jobinterval=5000 --disable-ui --nopreflight && cd ../../../"
    	fi > freebayes_02_FBlongrangerAlign_single_${CONT_DB}.${slurmID}.plan                
        
        echo "$(${LONGRANGER_PATH}/longranger mkref --version)" > freebayes_02_FBlongrangerAlign_single_${CONT_DB}.${slurmID}.version
        echo "$(${LONGRANGER_PATH}/longranger align --version)" >> freebayes_02_FBlongrangerAlign_single_${CONT_DB}.${slurmID}.version
	### 03_FBfreebayes
	elif [[ ${currentStep} -eq 3 ]]
    then
	    ### clean up plans 
        for x in $(ls freebayes_03_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
		REFNAME=$(basename ${CT_FREEBAYES_REFFASTA})
        
        if [[ ! -f "${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/ref/${REFNAME}" ]]
        then
    		(>&2 echo "ERROR - cannot find reference fasta file \"${REFNAME}\" in dir \"${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/ref\"")
        	exit 1
   		fi
   		
   		if [[ ! -d "${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/bams" ]]
        then
        	(>&2 echo "ERROR - cannot access directory ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/bams!")
        	exit 1
   		fi
        
        outdir="${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/freebayes/"
   		if [[ ! -d "${outdir}" ]]
        then
        	(>&2 echo "ERROR - cannot access directory ${outdir}!")
        	exit 1
   		fi
   	
   	
   		ref=${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/ref/refdata-${REFNAME%.fasta}/fasta/genome.fa
   		
   		if [[ ! -f "${ref}.fai" ]]
        then
        	(>&2 echo "ERROR - cannot reference fasta index file ${ref}.fai!")
        	exit 1
   		fi
   		
   		bam=${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/bams/10x_${PROJECT_ID}_longrangerAlign/outs/possorted_bam.bam
   		if [[ ! -f ${bam} ]]
   		then 
   			(>&2 echo "ERROR - cannot final duplicate marked bam file: ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/bams/${PROJECT_ID}_final10x.bam!")
        	exit 1
   		fi 
   		
   		echo "$(awk -v bam=${bam} -v ref=${ref} -v out=${outdir} '{print "freebayes --bam "bam" --region "$1":1-"$2" -f "ref" | bcftools view --no-version -Ou -o "out$1":1-"$2".bcf"}' ${ref}.fai)" > freebayes_03_FBfreebayes_block_${CONT_DB}.${slurmID}.plan

		echo "freebayes $(${CONDA_BASE_ENV} && freebayes --version && conda deactivate)" > freebayes_03_FBfreebayes_block_${CONT_DB}.${slurmID}.version
		echo "bcftools $(${CONDA_BASE_ENV} && bcftools --version | head -n1 | awk '{print $2}' && conda deactivate)" >> freebayes_03_FBfreebayes_block_${CONT_DB}.${slurmID}.version   
	### 04-FBconsensus
	elif [[ ${currentStep} -eq 4 ]]
    then
		### clean up plans 
        for x in $(ls freebayes_04_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        REFNAME=$(basename ${CT_FREEBAYES_REFFASTA})
        
		ref=${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/ref/refdata-${REFNAME%.fasta}/fasta/genome.fa
   		
   		if [[ ! -f "${ref}.fai" ]]
        then
        (>&2 echo "ERROR - cannot find reference fasta index file ${ref}.fai!")
        	exit 1
   		fi
   		
   		indir="${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/freebayes/"
   		if [[ ! -d "${indir}" ]]
        then
        	(>&2 echo "ERROR - cannot access directory ${indir}!")
        	exit 1
   		fi
   		
   		outdir="${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/"
        
    	# create list of bcf files, same order as in ref.fai
        echo "awk -v d=\"${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/freebayes/\" '{print d\$1\":1-\"\$2\".bcf\"}' ${ref}.fai > ${outdir}/${PROJECT_ID}_10x_concatList.txt" > freebayes_04_FBconsensus_single_${CONT_DB}.${slurmID}.plan
        echo "bcftools concat -Ou -f ${outdir}/${PROJECT_ID}_10x_concatList.txt | bcftools view -Ou -e'type=\"ref\"' | bcftools norm -Ob -f $ref -o ${outdir}/${PROJECT_ID}_10x.bcf" >> freebayes_04_FBconsensus_single_${CONT_DB}.${slurmID}.plan
        echo "bcftools index ${outdir}/${PROJECT_ID}_10x.bcf" >> freebayes_04_FBconsensus_single_${CONT_DB}.${slurmID}.plan
        echo "bcftools consensus -i'QUAL>1 && (GT=\"AA\" || GT=\"Aa\")' -Hla -f ${ref} ${outdir}/${PROJECT_ID}_10x.bcf > ${outdir}/${PROJECT_ID}_10x.fasta" >> freebayes_04_FBconsensus_single_${CONT_DB}.${slurmID}.plan
      
		echo "echo \"Num. bases affected \$(bcftools view -H -i 'QUAL>1 && (GT=\"AA\" || GT=\"Aa\")' -Ov ${outdir}/${PROJECT_ID}_10x.bcf | awk -F \"\\t\" '{print \$4\"\\t\"\$5}' | awk '{lenA=length(\$1); lenB=length(\$2); if (lenA < lenB ) {sum+=lenB-lenA} else if ( lenA > lenB ) { sum+=lenA-lenB } else {sum+=lenA}} END {print sum}')\" > ${outdir}/${PROJECT_ID}_10x.numvar " >> freebayes_04_FBconsensus_single_${CONT_DB}.${slurmID}.plan
		echo "bcftools view -i 'QUAL>1 && (GT=\"AA\" || GT=\"Aa\")' -Oz  ${outdir}/${PROJECT_ID}_10x.bcf > ${outdir}/${PROJECT_ID}_10x.changes.vcf.gz" >> freebayes_04_FBconsensus_single_${CONT_DB}.${slurmID}.plan
            
      	echo "bcftools $(${CONDA_BASE_ENV} && bcftools --version | head -n1 | awk '{print $2}' && conda deactivate)" > freebayes_04_FBconsensus_single_${CONT_DB}.${slurmID}.version
	### 05-FBstatistics
	elif [[ ${currentStep} -eq 5 ]]
    then
		### clean up plans 
        for x in $(ls freebayes_05_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
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
    	echo "${SUBMIT_SCRIPTS_PATH}/assemblyStats.sh ${configFile} 11" > freebayes_05_FBstatistics_single_${CONT_DB}.${slurmID}.plan
    	git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD > freebayes_05_FBstatistics_single_${CONT_DB}.${slurmID}.version
	else
        (>&2 echo "step ${currentStep} in CT_FREEBAYES_TYPE ${CT_FREEBAYES_TYPE} not supported")
        (>&2 echo "valid steps are: ${myTypes[${CT_FREEBAYES_TYPE}]}")
        exit 1            
    fi    	    		
else
    (>&2 echo "unknown CT_FREEBAYES_TYPE ${CT_FREEBAYES_TYPE}")
    (>&2 echo "supported types")
    x=0; while [ $x -lt ${#myTypes[*]} ]; do (>&2 echo "${myTypes[${x}]}"); done 
    exit 1
fi

exit 0