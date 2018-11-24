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


# Type: 0 - for 10x data 
# Type: 1 - for Illuimina PE data
myTypes=("1-FBprepareInput, 2-FBfastp, 3-FBbwa, 4-FBmarkDuplicates, 5-FBfreebayes, 6-FBbcftools", 
"1-FBprepareInput, 2-FBbwa, 3-FBmarkDuplicates, 4-FBfreebayes, 5-FBbcftools")
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
		
		echo "samtools $(${PACBIO_BASE_ENV} && samtools 2>&1 | grep Version | awk '{print $2}' && ${PACBIO_BASE_ENV_DEACT})" > freebayes_01_FBprepareInput_single_${CONT_DB}.${slurmID}.version
		echo "bwa $(${PACBIO_BASE_ENV} && bwa 2>&1 | grep Version | awk '{print $2}' && ${PACBIO_BASE_ENV_DEACT})" >> freebayes_01_FBprepareInput_single_${CONT_DB}.${slurmID}.version
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
	   		echo "fastp $(${PACBIO_BASE_ENV} && fastp 2>&1 | grep version | awk '{print $2}' && ${PACBIO_BASE_ENV_DEACT})" > freebayes_02_FBfastp_block_${CONT_DB}.${slurmID}.version	   		
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
			
	   		echo "bwa $(${PACBIO_BASE_ENV} && bwa 2>&1 | grep Version | awk '{print $2}' && ${PACBIO_BASE_ENV_DEACT})" > freebayes_03_FBbwa_block_${CONT_DB}.${slurmID}.version
        	echo "samtools $(${PACBIO_BASE_ENV} && samtools 2>&1 | grep Version | awk '{print $2}' && ${PACBIO_BASE_ENV_DEACT})" >> freebayes_03_FBbwa_block_${CONT_DB}.${slurmID}.version
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
   		   		   				
		for ib in ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/bams/*_bwa.bam
		do
			d=$(dirname ${ib})
			ob="${ib%_bwa.bam}_bwaMarkDups.bam"
			m="${ib%_bwa.bam}_bwaMarkDups.metrics"
			
			echo "picard MarkDuplicates I=${ib} O=${ob} M=${m} && samtools index -@ ${CT_FREEBAYES_SAMTOOLS_THREADS} ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/bams/${ob}" 				 
		done > freebayes_04_FBmarkDuplicates_block_${CONT_DB}.${slurmID}.plan
		
	   	echo "picard MarkDuplicates $(${PACBIO_BASE_ENV} && picard MarkDuplicates --version && ${PACBIO_BASE_ENV_DEACT})" > freebayes_04_FBmarkDuplicates_block_${CONT_DB}.${slurmID}.version
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
   		
   		### create list of input bam files 
   		numFiles=0
   		bamList=${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/bams/bamlist.txt
   		for ib in ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/bams/*_bwaMarkDups.bam
		do
			if [[ -f ${ib} ]]
			then
				echo "${ib}"	
				numFiles=$((${numFiles}+1))
			fi			
		done > ${bamList}
		
		if [[ ${numFiles} -eq 0 ]]
        then
    		(>&2 echo "ERROR - No duplicate marked bam files available!")
        	exit 1
   		fi
		
   		echo "$(awk -v bam=${bamList} -v ref=${ref} -v out=${outdir} '{print "freebayes --bam-list "bam" --region "$1":1-"$2" -f "ref" | bcftools view --no-version -Ob -o "out$1":1-"$2".bcf"}' ${ref}.fai)" > freebayes_05_FBfreebayes_block_${CONT_DB}.${slurmID}.plan

		echo "freebayes $(${PACBIO_BASE_ENV} && freebayes --version && ${PACBIO_BASE_ENV_DEACT})" > freebayes_05_FBfreebayes_block_${CONT_DB}.${slurmID}.version
		echo "bcftools $(${PACBIO_BASE_ENV} && bcftools --version | head -n1 | awk '{print $2}' && ${PACBIO_BASE_ENV_DEACT})" >> freebayes_05_FBfreebayes_block_${CONT_DB}.${slurmID}.version
   	### 6-FBbcftools
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
        echo "awk '{print \$1\":1-\"\$2\".bcf\"}' ${ref}.fai > ${outdir}/${PROJECT_ID}_10x_concatList.txt" > freebayes_06_FBconsensus_single_${CONT_DB}.${slurmID}.plan
        echo "bcftools concat -nf ${outdir}/${PROJECT_ID}_10x_concatList.txt | bcftools view -Ou -e'type="ref"' | bcftools norm -Ob -f $ref -o ${outdir}/${PROJECT_ID}_10x.bcf" >> freebayes_06_FBconsensus_single_${CONT_DB}.${slurmID}.plan
        echo "bcftools index ${outdir}/${PROJECT_ID}_10x.bcf" >> freebayes_06_FBconsensus_single_${CONT_DB}.${slurmID}.plan
        echo "bcftools consensus -i'QUAL>1 && (GT="AA" || GT="Aa")' -Hla -f ${ref} ${outdir}/${PROJECT_ID}_10x.bcf > ${outdir}/${PROJECT_ID}_10x.fasta" >> freebayes_06_FBconsensus_single_${CONT_DB}.${slurmID}.plan
      
      	echo "bcftools $(${PACBIO_BASE_ENV} && bcftools --version | head -n1 | awk '{print $2}' && ${PACBIO_BASE_ENV_DEACT})" > freebayes_06_FBconsensus_single_${CONT_DB}.${slurmID}.version  
                	   						
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