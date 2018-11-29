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

function setbwaOptions()
{
	CONTIG_BWA_OPT=""
	
	if [[ -z "${CT_HIC_BWA_THREADS}" ]]
	then 
		CT_HIC_BWA_THREADS=4	
	fi
	
	CONTIG_BWA_OPT="${CONTIG_BWA_OPT} -t ${CT_HIC_BWA_THREADS}"
	
	if [[ -n ${CT_HIC_BWA_VERBOSITY} ]]
	then 
		CONTIG_BWA_OPT="${CONTIG_BWA_OPT} -v ${CT_HIC_BWA_VERBOSITY}"
	fi
	
	if [[ -n ${CT_HIC_BWA_MISMATCHPENALTY} && ${CT_HIC_BWA_MISMATCHPENALTY} -gt 0 ]]
	then 
		CONTIG_BWA_OPT="${CONTIG_SAMTOOLS_OPT} -B ${CT_HIC_BWA_MISMATCHPENALTY}"
	fi
}


function setSamtoolsOptions()
{
	CONTIG_SAMTOOLS_OPT=""
	
	if [[ -z "${CT_HIC_SAMTOOLS_THREADS}" ]]
	then 
		CT_HIC_SAMTOOLS_THREADS=4	
	fi
	
	CONTIG_SAMTOOLS_OPT="${CONTIG_SAMTOOLS_OPT} -@ ${CT_HIC_SAMTOOLS_THREADS}"
	
	if [[ -n ${CT_HIC_SAMTOOLS_MEM} ]]
	then 
		CONTIG_SAMTOOLS_OPT="${CONTIG_SAMTOOLS_OPT} -m ${CT_HIC_SAMTOOLS_MEM}G"
	else
		CT_HIC_SAMTOOLS_MEM=2
		CONTIG_SAMTOOLS_OPT="${CONTIG_SAMTOOLS_OPT} -m ${CT_HIC_SAMTOOLS_MEM}G"
	fi		
}


# Type: 0 - Arima Mapping Pipeline (For QC) 
# Type: 1 - Phase Genomics Mapping Pipeline (For QC)
# Type: 2 - Aiden Lab Juicer Pipeline (For QC)
# Type: 3 - Salsa2 Pipeline (For Scaffolding)
# Type: 4 - 3d-dna Pipeline (For Scaffolding)
myTypes=("1-HICprepareInput, 2-HICbwa, 3-HICfilter, 4-HICmerge, 5-HICmarkduplicates, 6-HICstatistics", 
"1-HICprepareInput, 2-HICbwa, 3-HICfilter, 4-HICmatlock")
if [[ ${CT_HIC_TYPE} -eq 0 ]]
then 
    ### 1-HICprepareInput
    if [[ ${currentStep} -eq 1 ]]
    then
        ### clean up plans 
        for x in $(ls hic_01_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        if [[ ! -f "${CT_HIC_REFFASTA}" ]]
        then
        	(>&2 echo "ERROR - set CT_HIC_REFFASTA to reference fasta file")
        	exit 1
   		fi
   		
   		if [[ ! -d "${CT_HIC_READS}" ]]
        then
        	(>&2 echo "ERROR - set CT_HIC_READS to HiC read directory")
        	exit 1
   		fi
   		
   		echo "if [[ -d ${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID} ]]; then mv ${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID} ${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID}_$(date '+%Y-%m-%d_%H-%M-%S'); fi && mkdir ${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID}" > hic_01_HICprepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "mkdir -p ${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID}/reads" >> hic_01_HICprepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "mkdir -p ${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID}/bams" >> hic_01_HICprepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "mkdir -p ${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID}/ref" >> hic_01_HICprepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "mkdir -p ${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID}/config" >> hic_01_HICprepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "ln -s -r ${CT_HIC_REFFASTA} ${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID}/ref" >> hic_01_HICprepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "ln -s -r ${CT_HIC_REFFASTA} ${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID}/ref" >> hic_01_HICprepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "samtools faidx ${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID}/ref/$(basename ${CT_HIC_REFFASTA})" >> hic_01_HICprepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "bwa index ${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID}/ref/$(basename ${CT_HIC_REFFASTA})" >> hic_01_HICprepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "cp ${configFile} ${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID}/config/$(basename ${configFile%.sh})_$(date '+%Y-%m-%d_%H-%M-%S').sh" >> hic_01_HICprepareInput_single_${CONT_DB}.${slurmID}.plan
		
		echo "samtools $(${PACBIO_BASE_ENV} && samtools 2>&1 | grep Version | awk '{print $2}' && ${PACBIO_BASE_ENV_DEACT})" > hic_01_HICprepareInput_single_${CONT_DB}.${slurmID}.version
		echo "bwa $(${PACBIO_BASE_ENV} && bwa 2>&1 | grep Version | awk '{print $2}' && ${PACBIO_BASE_ENV_DEACT})" >> hic_01_HICprepareInput_single_${CONT_DB}.${slurmID}.version
	### 2-HICbwa 
    elif [[ ${currentStep} -eq 2 ]]
    then
        ### clean up plans 
        for x in $(ls hic_02_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
		
		if [[ ! -d "${CT_HIC_READS}" ]]
        then
        	(>&2 echo "ERROR - set CT_HIC_READS to directory that contain the PROJECT_ID*.fastq.qz read files")
        	exit 1
   		fi   		   				
        
        if [[ ! -d "${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID}/reads" ]]
        then
        	(>&2 echo "ERROR - cannot access directory ${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID}/reads!")
        	exit 1
   		fi
   		
   		ref=${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID}/ref/$(basename ${CT_HIC_REFFASTA})
   		
   		if [[ ! -f "${ref}" ]]
        then
        (>&2 echo "ERROR - cannot reference fasta file ${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID}/ref/$(basename ${CT_HIC_REFFASTA})!")
        	exit 1
   		fi
   		
   		### link HiC reads into current reads sub directory
   		for x in ${CT_HIC_READS}/${PROJECT_ID}_*_*_R[12].fastq.gz
		do
			if [[ -f ${x} ]]
			then	
				ln -s -r -f ${x} ${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID}/reads 
			fi
		done
   		   				
   		numR1Files=0
		for x in ${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID}/reads/${PROJECT_ID}_*_*_R1.fastq.gz
		do
			if [[ -f ${x} ]]
			then	
				numR1Files=$((${numR1Files}+1))	
			fi
		done
		
		if [[ ${numR1Files} -eq 0 ]]
        then
        	(>&2 echo "ERROR - cannot read HiC R1 files with following pattern: ${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID}/reads/${PROJECT_ID}_*_*_R1.fastq.gz")
        	exit 1
   		fi
   		
   		numR2Files=0
		for x in ${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID}/reads/${PROJECT_ID}_*_*_R2.fastq.gz
		do
			if [[ -f ${x} ]]
			then	
				numR2Files=$((${numR2Files}+1))	
			fi
		done
		
		if [[ ${numR2Files} -eq 0 ]]
        then
        	(>&2 echo "ERROR - cannot read HiC R2 files with following pattern: ${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID}/reads/${PROJECT_ID}_*_*_R2.fastq.gz")
        	exit 1
   		fi
   		
   		if [[ ${numR1Files} -ne ${numR2Files} ]]
        then
        	(>&2 echo "ERROR - HiC R1 files ${numR1Files} does not match R2 files ${numR2Files}")
        	exit 1
   		fi
   		
   		setbwaOptions
   		
		for r1 in ${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID}/reads/${PROJECT_ID}_*_*_R1.fastq.gz
		do
			id=$(dirname ${r1})
			f1=$(basename ${r1})
			f2=$(echo "${f1}" | sed -e "s:_R1.fastq.gz:_R2.fastq.gz:")
			o="${f1%_R1.fastq.gz}"											
			
			echo "bwa mem${CONTIG_BWA_OPT} -R \"@RG\tID:${o}\tSM:${PROJECT_ID}_HIC\tLB:${PROJECT_ID}_HIC\tPL:ILLUMINA:\tPU:none\" ${ref} ${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID}/reads/${f1} | samtools view -Sb - > ${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID}/bams/${o}_bwa_1.bam"
			echo "bwa mem${CONTIG_BWA_OPT} -R \"@RG\tID:${o}\tSM:${PROJECT_ID}_HIC\tLB:${PROJECT_ID}_HIC\tPL:ILLUMINA:\tPU:none\" ${ref} ${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID}/reads/${f2} | samtools view -Sb - > ${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID}/bams/${o}_bwa_2.bam" 				 
		done > hic_02_HICbwa_block_${CONT_DB}.${slurmID}.plan
		
   		echo "bwa $(${PACBIO_BASE_ENV} && bwa 2>&1 | grep Version | awk '{print $2}' && ${PACBIO_BASE_ENV_DEACT})" > hic_02_HICbwa_block_${CONT_DB}.${slurmID}.version
   		echo "samtools $(${PACBIO_BASE_ENV} && samtools 2>&1 | grep Version | awk '{print $2}' && ${PACBIO_BASE_ENV_DEACT})" >> hic_02_HICbwa_block_${CONT_DB}.${slurmID}.version
	### 3-HICfilter
    elif [[ ${currentStep} -eq 3 ]]
    then
        ### clean up plans 
        for x in $(ls hic_03_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done

		if [[ ! -d "${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID}/bams" ]]
        then
        	(>&2 echo "ERROR - cannot access directory ${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID}/bams!")
        	exit 1
   		fi
   		   		   				
		for b1 in ${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID}/bams/*_bwa_1.bam
		do
			d=$(dirname ${b1})
			b2="${b1%_1.bam}_2.bam"
			f1=$(basename ${b1%_bwa_1.bam})_bwaFilt_1.bam
			f2="${f1%_1.bam}_2.bam"			
			
			echo "samtools view -h ${b1} | perl ${MARVEL_PATH}/scripts/filter_five_end.pl | samtools view -Sb - > ${d}/${f1}"
			echo "samtools view -h ${b2} | perl ${MARVEL_PATH}/scripts/filter_five_end.pl | samtools view -Sb - > ${d}/${f2}" 				 
		done > hic_03_HICfilter_block_${CONT_DB}.${slurmID}.plan
					
		echo "samtools $(${PACBIO_BASE_ENV} && samtools 2>&1 | grep Version | awk '{print $2}' && ${PACBIO_BASE_ENV_DEACT})" > hic_03_HICfilter_block_${CONT_DB}.${slurmID}.version	   	
	### 4-HICmerge
    elif [[ ${currentStep} -eq 4 ]]
    then
        ### clean up plans 
        for x in $(ls hic_04_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        if [[ ! -d "${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID}/bams" ]]
        then
        	(>&2 echo "ERROR - cannot access directory ${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID}/bams!")
        	exit 1
   		fi
   		
   		ref=${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID}/ref/$(basename ${CT_HIC_REFFASTA})
   		
   		if [[ ! -f ${ref}.fai ]]
   		then  
   		 	(>&2 echo "ERROR - cannot access reference fasta index ${ref}.fai!")
        	exit 1
		fi
		 
		if [[ -z ${CT_HIC_MINMAPQV} ]]
		then
			CT_HIC_MINMAPQV=10	
		fi
		 
		for b1 in ${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID}/bams/*_bwaFilt_1.bam
		do
			b2="${b1%_1.bam}_2.bam"
			o="${b1%_1.bam}.bam"			
			
			echo "perl ${MARVEL_PATH}/scripts/two_read_bam_combiner.pl ${b1} ${b2} $(which samtools) ${CT_HIC_MINMAPQV} | samtools view -bS -t ${ref}.fai - | samtools sort -o ${o} -"			 				 
			done > hic_04_HICmerge_single_${CONT_DB}.${slurmID}.plan
		   
		echo "samtools $(${PACBIO_BASE_ENV} && samtools 2>&1 | grep Version | awk '{print $2}' && ${PACBIO_BASE_ENV_DEACT})" > hic_03_HICfilter_single_${CONT_DB}.${slurmID}.version		
   	### 5-HICmarkduplicates
    elif [[ ${currentStep} -eq 5 ]]
    then
        ### clean up plans 
        for x in $(ls hic_05_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
		if [[ ! -d "${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID}/bams" ]]
        then
        	(>&2 echo "ERROR - cannot access directory ${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID}/bams!")
        	exit 1
   		fi
   		   		
   		## if multiple bam files are available (e.g. different Lanes) then merge files prior to markduplicates
   		files=$(ls ${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID}/bams/*_bwaFilt.bam)
   		
   		echo "picard MarkDuplicates $(${PACBIO_BASE_ENV} && picard MarkDuplicates --version && ${PACBIO_BASE_ENV_DEACT})" > hic_05_HICmarkduplicates_single_${CONT_DB}.${slurmID}.version
   		
   		if [[ $(echo $files | wc -w) -eq 1 ]]
   		then
   			ob="${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID}/bams/${PROJECT_ID}_finalHiC.bam"
			m="${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID}/bams/${PROJECT_ID}_finalHiC.metrics"
		   echo "picard MarkDuplicates I=${files} O=${ob} M=${m} && samtools index -@ ${CT_HIC_SAMTOOLS_THREADS} ${ob} && ln -s -f -r ${ob} ${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID} && ln -s -f -r ${ob}.bai ${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID}"
   		elif [[ $(echo $files | wc -w) -gt 1 ]]
   		then
   			mrg=${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID}/bams/${PROJECT_ID}_mergedHiC.bam
   			o=${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID}/bams/${PROJECT_ID}_finalHiC.bam
   			m=${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID}/bams/${PROJECT_ID}_finalHiC.metrics
   			i=$(echo -e ${files} | sed -e "s:${CT_HIC_OUTDIR}:I=${CT_HIC_OUTDIR}:g")
   			echo "picard MergeSamFiles ${i} OUTPUT=${mrg} USE_THREADING=TRUE ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT && picard MarkDuplicates I=${mrg} O=${o} M=${m} && samtools index -@ ${CT_HIC_SAMTOOLS_THREADS} ${o} && ln -s -f -r ${o} ${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID} && ln -s -f -r ${ob}.bai ${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID}"
   			echo "picard MergeSamFiles $(${PACBIO_BASE_ENV} && picard MarkDuplicates --version && ${PACBIO_BASE_ENV_DEACT})" >> hic_05_HICmarkduplicates_single_${CONT_DB}.${slurmID}.version	
   		else
   	 		(>&2 echo "ERROR - cannot find file with following pattern: ${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID}/bams/*_bwaFilt.bam!")
        	exit 1
   	 	fi > hic_05_HICmarkduplicates_single_${CONT_DB}.${slurmID}.plan
   	### 6-HICstatistics
    elif [[ ${currentStep} -eq 6 ]]
    then
        ### clean up plans 
        for x in $(ls hic_06_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
		if [[ ! -f "${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID}/bams/${PROJECT_ID}_finalHiC.bam" ]]
        then
    	(>&2 echo "ERROR - cannot access final duplicate marked bam file ${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID}/bams/${PROJECT_ID}_finalHiC.bam!")
        	exit 1
   		fi
   		   		  
   		echo "perl ${MARVEL_PATH}/scripts/get_stats.pl ${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID}/bams/${PROJECT_ID}_finalHiC.bam > ${CT_HIC_OUTDIR}/hic_${CT_HIC_RUNID}/bams/${PROJECT_ID}_finalHiC.stats" > hic_06_HICstatistics _single_${CONT_DB}.${slurmID}.plan
    else
        (>&2 echo "step ${currentStep} in CT_HIC_TYPE ${CT_HIC_TYPE} not supported")
        (>&2 echo "valid steps are: ${myTypes[${CT_HIC_TYPE}]}")
        exit 1            
    fi    		
else
    (>&2 echo "unknown CT_HIC_TYPE ${CT_HIC_TYPE}")
    (>&2 echo "supported types")
    x=0; while [ $x -lt ${#myTypes[*]} ]; do (>&2 echo "${myTypes[${x}]}"); done 
    exit 1
fi

exit 0