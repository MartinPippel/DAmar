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
myTypes=("1-WHprepareInput, 2-WHminimap2, 3-")
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
   		
   		echo "if [[ -d ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID} ]]; then mv ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID} ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}_$(date '+%Y-%m-%d_%H-%M-%S'); fi && mkdir ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}" > whatshap_01_WHprepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "mkdir -p ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/pb_reads" >> whatshap_01_WHprepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "mkdir -p ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/10x_reads" >> whatshap_01_WHprepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "mkdir -p ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/pb_bams" >> whatshap_01_WHprepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "mkdir -p ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/10x_bams" >> whatshap_01_WHprepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "mkdir -p ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/ref" >> whatshap_01_WHprepareInput_single_${CONT_DB}.${slurmID}.plan		
		echo "ln -s -r ${CT_WHATSHAP_REFFASTA} ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/ref" >> whatshap_01_WHprepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "samtools faidx ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/ref/$(basename ${CT_WHATSHAP_REFFASTA})" >> whatshap_01_WHprepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "bwa index ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/ref/$(basename ${CT_WHATSHAP_REFFASTA})" >> whatshap_01_WHprepareInput_single_${CONT_DB}.${slurmID}.plan
		
		echo "samtools $(${PACBIO_BASE_ENV} && samtools 2>&1 | grep Version | awk '{print $2}' && ${PACBIO_BASE_ENV_DEACT})" > whatshap_01_WHprepareInput_single_${CONT_DB}.${slurmID}.version
		echo "bwa $(${PACBIO_BASE_ENV} && bwa 2>&1 | grep Version | awk '{print $2}' && ${PACBIO_BASE_ENV_DEACT})" >> whatshap_01_WHprepareInput_single_${CONT_DB}.${slurmID}.version
	### 2-WHminimap2 
    elif [[ ${currentStep} -eq 2 ]]
    then
        ### clean up plans 
        for x in $(ls whatshap_02_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
		
		if [[ ! -d "${CT_WHATSHAP_READS}" ]]
        then
        	(>&2 echo "ERROR - set CT_WHATSHAP_READS to directory that contain the PROJECT_ID*.fastq.qz read files")
        	exit 1
   		fi
   		   				
		# we need to trim off barcode and adapter for 10x 
		if [[ "x${CT_WHATSHAP_READSTYPE}" == "x10x" ]]
		then
	   		numR1Files=0
			for x in ${CT_WHATSHAP_READS}/${PROJECT_ID}_S*_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz
			do
				if [[ -f ${x} ]]
				then	
					numR1Files=$((${numR1Files}+1))	
				fi
			done
			
			if [[ ${numR1Files} -eq 0 ]]
	        then
	        	(>&2 echo "ERROR - cannot read 10x R1 files with following pattern: ${CT_WHATSHAP_READS}/${PROJECT_ID}_S*_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz")
	        	exit 1
	   		fi
	   		
	   		numR2Files=0
			for x in ${CT_WHATSHAP_READS}/${PROJECT_ID}_S*_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz
			do
				if [[ -f ${x} ]]
				then	
					numR2Files=$((${numR2Files}+1))	
				fi
			done
			
			if [[ ${numR2Files} -eq 0 ]]
	        then
	        	(>&2 echo "ERROR - cannot read 10x R2 files with following pattern: ${CT_WHATSHAP_READS}/${PROJECT_ID}_S*_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz")
	        	exit 1
	   		fi
	   		
	   		if [[ ${numR1Files} -ne ${numR2Files} ]]
	        then
	        	(>&2 echo "ERROR - 10x R1 files ${numR1Files} does not match R2 files ${numR2Files}")
	        	exit 1
	   		fi
	   		
	   		setFastpOptions
	   		
			for r1 in ${CT_WHATSHAP_READS}/${PROJECT_ID}_S[0-9]_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz
			do
				id=$(dirname ${r1})
				f1=$(basename ${r1})
				f2=$(echo "${f1}" | sed -e "s:_R1_:_R2_:")
				o="${f1%_R1_???.fastq.gz}"											
				
				echo "fastp -i ${id}/${f1} -I ${id}/${f2} -f 23 -G -Q -j ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/reads/${o}.json -h ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/reads/${o}.html -w ${CT_WHATSHAP_FASTP_THREADS} -o ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/reads/${f1} -O ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/reads/${f2}"				 
			done > whatshap_02_FBfastp_block_${CONT_DB}.${slurmID}.plan
	   		echo "fastp $(${PACBIO_BASE_ENV} && fastp 2>&1 | grep version | awk '{print $2}' && ${PACBIO_BASE_ENV_DEACT})" > whatshap_02_FBfastp_block_${CONT_DB}.${slurmID}.version	   		
		else		
			(>&2 echo "ERROR - set CT_WHATSHAP_READSTYPE ${CT_WHATSHAP_READSTYPE} not supported yet")
	        exit 1			
		fi
	### 3-FBbwa
    elif [[ ${currentStep} -eq 3 ]]
    then
    	### clean up plans 
        for x in $(ls whatshap_03_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        if [[ ! -d "${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/reads" ]]
        then
        	(>&2 echo "ERROR - cannot access directory ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/reads!")
        	exit 1
   		fi
   		
   		ref=${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/ref/$(basename ${CT_WHATSHAP_REFFASTA})
   		
   		if [[ ! -f "${ref}" ]]
        then
        (>&2 echo "ERROR - cannot reference fasta file ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/ref/$(basename ${CT_WHATSHAP_REFFASTA})!")
        	exit 1
   		fi
   		   				
		# we need to trim off barcode and adapter for 10x 
		if [[ "x${CT_WHATSHAP_READSTYPE}" == "x10x" ]]
		then
	   		numR1Files=0
			for x in ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/reads/${PROJECT_ID}_S*_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz
			do
				if [[ -f ${x} ]]
				then	
					numR1Files=$((${numR1Files}+1))	
				fi
			done
			
			if [[ ${numR1Files} -eq 0 ]]
	        then
	        	(>&2 echo "ERROR - cannot read 10x R1 files with following pattern: ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/reads/${PROJECT_ID}_S*_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz")
	        	exit 1
	   		fi
	   		
	   		numR2Files=0
			for x in ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/reads/${PROJECT_ID}_S*_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz
			do
				if [[ -f ${x} ]]
				then	
					numR2Files=$((${numR2Files}+1))	
				fi
			done
			
			if [[ ${numR2Files} -eq 0 ]]
	        then
	        	(>&2 echo "ERROR - cannot read 10x R2 files with following pattern: ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/reads/${PROJECT_ID}_S*_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz")
	        	exit 1
	   		fi
	   		
	   		if [[ ${numR1Files} -ne ${numR2Files} ]]
	        then
	        	(>&2 echo "ERROR - 10x R1 files ${numR1Files} does not match R2 files ${numR2Files}")
	        	exit 1
	   		fi
	   		
	   		setbwaOptions
	   		setSamtoolsOptions
	   		
			for r1 in ${CT_WHATSHAP_READS}/${PROJECT_ID}_S[0-9]_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz
			do
				id=$(dirname ${r1})
				f1=$(basename ${r1})
				f2=$(echo "${f1}" | sed -e "s:_R1_:_R2_:")
				o="${f1%_R1_???.fastq.gz}"											
				
				echo "bwa mem${CONTIG_BWA_OPT} -R \"@RG\tID:${PROJECT_ID}\tSM:10X\" ${ref} ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/reads/${f1} ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/reads/${f2} | samtools sort${CONTIG_SAMTOOLS_OPT} -O BAM -o ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/bams/${o}_bwa.bam && samtools index -@ ${CT_WHATSHAP_SAMTOOLS_THREADS} ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/bams/${o}_bwa.bam" 				 
			done > whatshap_03_FBbwa_block_${CONT_DB}.${slurmID}.plan
			
	   		echo "bwa $(${PACBIO_BASE_ENV} && bwa 2>&1 | grep Version | awk '{print $2}' && ${PACBIO_BASE_ENV_DEACT})" > whatshap_03_FBbwa_block_${CONT_DB}.${slurmID}.version
        	echo "samtools $(${PACBIO_BASE_ENV} && samtools 2>&1 | grep Version | awk '{print $2}' && ${PACBIO_BASE_ENV_DEACT})" >> whatshap_03_FBbwa_block_${CONT_DB}.${slurmID}.version
		else		
			(>&2 echo "ERROR - set CT_WHATSHAP_READSTYPE ${CT_WHATSHAP_READSTYPE} not supported yet")
	        exit 1			
		fi        	
    ### 4-FBmarkDuplicates
    elif [[ ${currentStep} -eq 4 ]]
    then
        ### clean up plans 
        for x in $(ls whatshap_04_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done               

		if [[ ! -d "${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/bams" ]]
        then
        	(>&2 echo "ERROR - cannot access directory ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/bams!")
        	exit 1
   		fi
   		
   		## if multiple bam files are available (e.g. different Lanes) then merge files prior to markduplicates
   		files=$(ls ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/bams/*_bwa.bam)
   		
   		echo "picard MarkDuplicates $(${PACBIO_BASE_ENV} && picard MarkDuplicates --version && ${PACBIO_BASE_ENV_DEACT})" > whatshap_04_FBmarkDuplicates_block_${CONT_DB}.${slurmID}.version
   		if [[ $(echo $files | wc -w) -eq 1 ]]
   		then
   			ob="${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/bams/${PROJECT_ID}_final10x.bam"
			m="${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/bams/${PROJECT_ID}_final10x.metrics"
   			
   			echo "picard MarkDuplicates I=${files} O=${ob} M=${m} && samtools index -@ ${CT_WHATSHAP_SAMTOOLS_THREADS} ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/bams/${ob}" 				 
   		elif [[ $(echo $files | wc -w) -gt 1 ]]
   		then
   			mrg=${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/bams/${PROJECT_ID}_merged10x.bam
   			o=${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/bams/${PROJECT_ID}_final10x.bam
   			m=${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/bams/${PROJECT_ID}_final10x.metrics
   			i=$(echo -e ${files} | sed -e "s:${CT_WHATSHAP_OUTDIR}:I=${CT_WHATSHAP_OUTDIR}:g")
   			echo "picard MergeSamFiles ${i} OUTPUT=${mrg} USE_THREADING=TRUE ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT && picard MarkDuplicates I=${mrg} O=${o} M=${m} && samtools index -@ ${CT_WHATSHAP_SAMTOOLS_THREADS} ${o}"
   			echo "picard MergeSamFiles $(${PACBIO_BASE_ENV} && picard MarkDuplicates --version && ${PACBIO_BASE_ENV_DEACT})" >> whatshap_04_FBmarkDuplicates_block_${CONT_DB}.${slurmID}.version	   			
		else
   	 		(>&2 echo "ERROR - cannot find file with following pattern: ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/bams/*_bwa.bam")
        	exit 1
   	 	fi > whatshap_04_FBmarkDuplicates_block_${CONT_DB}.${slurmID}.plan   			   		   		   						   	
	### 5-FBwhatshap
    elif [[ ${currentStep} -eq 5 ]]
    then
        ### clean up plans 
        for x in $(ls whatshap_05_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        if [[ ! -d "${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/bams" ]]
        then
        	(>&2 echo "ERROR - cannot access directory ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/bams!")
        	exit 1
   		fi
   		
   		outdir="${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/whatshap/"
   		if [[ ! -d "${outdir}" ]]
        then
        	(>&2 echo "ERROR - cannot access directory ${outdir}!")
        	exit 1
   		fi   		
        
        ref=${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/ref/$(basename ${CT_WHATSHAP_REFFASTA})
   		
   		if [[ ! -f "${ref}.fai" ]]
        then
        	(>&2 echo "ERROR - cannot reference fasta index file ${ref}.fai!")
        	exit 1
   		fi
   		
   		bam=${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/bams/${PROJECT_ID}_final10x.bam
   		if [[ ! -f ${bam} ]]
   		then 
   			(>&2 echo "ERROR - cannot final duplicate marked bam file: ${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/bams/${PROJECT_ID}_final10x.bam!")
        	exit 1
   		fi 
   		
   		echo "$(awk -v bam=${bam} -v ref=${ref} -v out=${outdir} '{print "whatshap --bam "bam" --region "$1":1-"$2" -f "ref" | bcftools view --no-version -Ob -o "out$1":1-"$2".bcf"}' ${ref}.fai)" > whatshap_05_FBwhatshap_block_${CONT_DB}.${slurmID}.plan

		echo "whatshap $(${PACBIO_BASE_ENV} && whatshap --version && ${PACBIO_BASE_ENV_DEACT})" > whatshap_05_FBwhatshap_block_${CONT_DB}.${slurmID}.version
		echo "bcftools $(${PACBIO_BASE_ENV} && bcftools --version | head -n1 | awk '{print $2}' && ${PACBIO_BASE_ENV_DEACT})" >> whatshap_05_FBwhatshap_block_${CONT_DB}.${slurmID}.version
   	### 6-FBbcftools
    elif [[ ${currentStep} -eq 6 ]]
    then
        ### clean up plans 
        for x in $(ls whatshap_06_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
		ref=${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/ref/$(basename ${CT_WHATSHAP_REFFASTA})
   		
   		if [[ ! -f "${ref}.fai" ]]
        then
        	(>&2 echo "ERROR - cannot reference fasta index file ${ref}.fai!")
        	exit 1
   		fi
   		
   		indir="${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/whatshap/"
   		if [[ ! -d "${indir}" ]]
        then
        	(>&2 echo "ERROR - cannot access directory ${indir}!")
        	exit 1
   		fi
   		
   		outdir="${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/"
        
    	# create list of bcf files, same order as in ref.fai
        echo "awk -v d=\"${CT_WHATSHAP_OUTDIR}/whatshap_${CT_WHATSHAP_RUNID}/whatshap/\" '{print d\$1\":1-\"\$2\".bcf\"}' ${ref}.fai > ${outdir}/${PROJECT_ID}_10x_concatList.txt" > whatshap_06_FBconsensus_single_${CONT_DB}.${slurmID}.plan
        echo "bcftools concat -nf ${outdir}/${PROJECT_ID}_10x_concatList.txt | bcftools view -Ou -e'type=\"ref\"' | bcftools norm -Ob -f $ref -o ${outdir}/${PROJECT_ID}_10x.bcf" >> whatshap_06_FBconsensus_single_${CONT_DB}.${slurmID}.plan
        echo "bcftools index ${outdir}/${PROJECT_ID}_10x.bcf" >> whatshap_06_FBconsensus_single_${CONT_DB}.${slurmID}.plan
        echo "bcftools consensus -i'QUAL>1 && (GT=\"AA\" || GT=\"Aa\")' -Hla -f ${ref} ${outdir}/${PROJECT_ID}_10x.bcf > ${outdir}/${PROJECT_ID}_10x.fasta" >> whatshap_06_FBconsensus_single_${CONT_DB}.${slurmID}.plan
      
      	echo "bcftools $(${PACBIO_BASE_ENV} && bcftools --version | head -n1 | awk '{print $2}' && ${PACBIO_BASE_ENV_DEACT})" > whatshap_06_FBconsensus_single_${CONT_DB}.${slurmID}.version  
                	   						
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