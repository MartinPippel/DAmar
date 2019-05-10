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
	
	if [[ -z "${SC_HIC_BWA_THREADS}" ]]
	then 
		SC_HIC_BWA_THREADS=4	
	fi
	
	CONTIG_BWA_OPT="${CONTIG_BWA_OPT} -t ${SC_HIC_BWA_THREADS}"
	
	if [[ -n ${SC_HIC_BWA_VERBOSITY} ]]
	then 
		CONTIG_BWA_OPT="${CONTIG_BWA_OPT} -v ${SC_HIC_BWA_VERBOSITY}"
	fi
	
	if [[ -n ${SC_HIC_BWA_MISMATCHPENALTY} && ${SC_HIC_BWA_MISMATCHPENALTY} -gt 0 ]]
	then 
		CONTIG_BWA_OPT="${CONTIG_BWA_OPT} -B ${SC_HIC_BWA_MISMATCHPENALTY}"
	fi
}


function setSamtoolsOptions()
{
	CONTIG_SAMTOOLS_OPT=""
	
	if [[ -z "${SC_HIC_SAMTOOLS_THREADS}" ]]
	then 
		SC_HIC_SAMTOOLS_THREADS=4	
	fi
	
	CONTIG_SAMTOOLS_OPT="${CONTIG_SAMTOOLS_OPT} -@ ${SC_HIC_SAMTOOLS_THREADS}"
	
	if [[ -n ${SC_HIC_SAMTOOLS_MEM} ]]
	then 
		CONTIG_SAMTOOLS_OPT="${CONTIG_SAMTOOLS_OPT} -m ${SC_HIC_SAMTOOLS_MEM}G"
	else
		SC_HIC_SAMTOOLS_MEM=2
		CONTIG_SAMTOOLS_OPT="${CONTIG_SAMTOOLS_OPT} -m ${SC_HIC_SAMTOOLS_MEM}G"
	fi		
}

function setPicardOptions()
{
	CONTIG_PICARD_OPT=""
	
	if [[ -z ${SC_HIC_PICARD_XMX} ]]
	then 
		SC_HIC_PICARD_XMX=8	
	fi
	
	if [[ -z ${SC_HIC_PICARD_XMS} ]]
	then 
		SC_HIC_PICARD_XMS=8	
	fi
	
	CONTIG_PICARD_OPT="-Xmx${SC_HIC_PICARD_XMX}G -Xms${SC_HIC_PICARD_XMS}G -XX:-UseGCOverheadLimit"	
}

function setJuicerOptions()
{
	SC_HIC_JUICER_OPT=""
	
	if [[ $# -eq 0 ]]
	then 
		ignoreStage=0	
	else
		ignoreStage=1
	fi
	
	# set juicer stage 
	if [[ -n ${SC_HIC_JUICER_STAGE} && ${ignoreStage} -eq 0 ]]
	then 
		
		if [[ ! "x${SC_HIC_JUICER_STAGE}" == "xmerge" && ! "x${SC_HIC_JUICER_STAGE}" == "xdedup" && ! "x${SC_HIC_JUICER_STAGE}" == "xfinal" && ! "x${SC_HIC_JUICER_STAGE}" == "xpostproc" && ! "x${SC_HIC_JUICER_STAGE}" == "xearly" ]]
		then
			(>&2 echo "[ERRROR] Unsupported juicer stage ${SC_HIC_JUICER_STAGE}. Can be: [merge, dedup, final, postproc, early]")
			exit 1	
		fi	
			
		SC_HIC_JUICER_OPT="${SC_HIC_JUICER_OPT} -S ${SC_HIC_JUICER_STAGE}"
	fi
	
	# set submission queue - short queue
	if [[ -n ${SC_HIC_JUICER_SHORTQUEUE} ]]
	then
		SC_HIC_JUICER_OPT="${SC_HIC_JUICER_OPT} -q ${SC_HIC_JUICER_SHORTQUEUE}"
	fi
	
	# set submission queue - short queue time limit, definition: time limit for queue, i.e. -W 12:00 is 12 hours, default (1200)
	if [[ -n ${SC_HIC_JUICER_SHORTQUEUETLIMIIT} ]]
	then
		SC_HIC_JUICER_OPT="${SC_HIC_JUICER_OPT} -Q ${SC_HIC_JUICER_SHORTQUEUETLIMIIT}"
	fi
	 
	# set submission queue - long queue
	if [[ -n ${SC_HIC_JUICER_LONGQUEUE} ]]
	then
		SC_HIC_JUICER_OPT="${SC_HIC_JUICER_OPT} -l ${SC_HIC_JUICER_LONGQUEUE}"
	fi
	
	# set submission queue - short queue time limit, definition: time limit for long queue, i.e. -W 168:00 is one week (default 3600)
	if [[ -n ${SC_HIC_JUICER_LONGQUEUETLIMIIT} ]]
	then
		SC_HIC_JUICER_OPT="${SC_HIC_JUICER_OPT} -L ${SC_HIC_JUICER_LONGQUEUETLIMIIT}"
	fi
	
	# chunk size definition: number of lines in split files, must be multiple of 4 (default 90000000, which equals 22.5 million reads)
	if [[ -n ${SC_HIC_JUICER_CHUNKSIZE} && ${SC_HIC_JUICER_CHUNKSIZE} -gt 0 ]]
	then
		SC_HIC_JUICER_OPT="${SC_HIC_JUICER_OPT} -C ${SC_HIC_JUICER_CHUNKSIZE}"
	fi
}

function setThreeDDNAOptions()
{
	THREEDDNA_OPT=""
	
	if [[ -n ${SC_HIC_3DDNA_MODE} ]]
	then 
		
		if [[ ! "x${SC_HIC_3DDNA_MODE}" == "xhaploid" && ! "x${SC_HIC_3DDNA_MODE}" == "xdiploid" ]]
		then
			(>&2 echo "[ERRROR] Unsupported 3d-dna mode ${SC_HIC_3DDNA_MODE}. Can be: [haploid, diploid]")
			exit 1	
		fi
		
		THREEDDNA_OPT="${THREEDDNA_OPT} -m ${SC_HIC_3DDNA_MODE}"
	fi	
	
	if [[ -n ${SC_HIC_3DDNA_MINCONTIGLEN} ]]
	then
		THREEDDNA_OPT="${THREEDDNA_OPT} -i ${SC_HIC_3DDNA_MINCONTIGLEN}"
	fi
	
	if [[ -n ${SC_HIC_3DDNA_ROUNDS} && ${SC_HIC_3DDNA_ROUNDS} -gt 0 ]]
	then
		THREEDDNA_OPT="${THREEDDNA_OPT} -r ${SC_HIC_3DDNA_ROUNDS}"
	fi
	
	if [[ -n ${SC_HIC_3DDNA_STAGE} ]]
	then
	if [[ ! "x${SC_HIC_3DDNA_STAGE}" == "xpolish" && ! "x${SC_HIC_3DDNA_STAGE}" == "xsplit" && ! "x${SC_HIC_3DDNA_STAGE}" == "xseal" && ! "x${SC_HIC_3DDNA_STAGE}" == "xmerge" && ! "x${SC_HIC_3DDNA_STAGE}" == "xfinalize" ]]
		then
			(>&2 echo "[ERRROR] Unsupported 3d-dna stage ${SC_HIC_3DDNA_STAGE}. can be: [polish, split, seal, merge, finalize]")
			exit 1	
		fi
		THREEDDNA_OPT="${THREEDDNA_OPT} -s ${SC_HIC_3DDNA_STAGE}"
	fi
	
	if [[ -n ${SC_HIC_3DDNA_MAPQV} ]]
	then
		THREEDDNA_OPT="${THREEDDNA_OPT} -q ${SC_HIC_3DDNA_MAPQV}"
	fi
}

function set3DDNAVisualizeOptions  
{
	THREEDDNA_VISUALIZE_OPT=""

	if [[ -n ${SC_HIC_3DDNAVISUALIZE_MAPQV} && ${SC_HIC_3DDNAVISUALIZE_MAPQV} -ge 0 && ${SC_HIC_3DDNAVISUALIZE_MAPQV} -le 100 ]]
	then 
		THREEDDNA_VISUALIZE_OPT="${THREEDDNA_VISUALIZE_OPT} -q ${SC_HIC_3DDNAVISUALIZE_MAPQV}"
	fi
	
	if [[ -n ${SC_HIC_3DDNAVISUALIZE_MNDPATH} ]]
	then 
		THREEDDNA_VISUALIZE_OPT="${THREEDDNA_VISUALIZE_OPT} -m ${SC_HIC_3DDNAVISUALIZE_MNDPATH}"
	fi
	
	if [[ -n ${SC_HIC_3DDNAVISUALIZE_SKIPNORM} && ${SC_HIC_3DDNAVISUALIZE_SKIPNORM} -gt 0 ]]
	then 
		THREEDDNA_VISUALIZE_OPT="${THREEDDNA_VISUALIZE_OPT} -n"
	fi
	
	if [[ -n ${SC_HIC_3DDNAVISUALIZE_RESOLUTION} ]]
	then 
		THREEDDNA_VISUALIZE_OPT="${THREEDDNA_VISUALIZE_OPT} -r ${SC_HIC_3DDNAVISUALIZE_RESOLUTION}"
	fi
	
	if [[ -n ${SC_HIC_3DDNAVISUALIZE_CLEANUP} && ${SC_HIC_3DDNAVISUALIZE_CLEANUP} -gt 0 ]]
	then 
		THREEDDNA_VISUALIZE_OPT="${THREEDDNA_VISUALIZE_OPT} -c"
	fi
	
	if [[ -n ${SC_HIC_3DDNAVISUALIZE_IGNOREMAPQV} && ${SC_HIC_3DDNAVISUALIZE_IGNOREMAPQV} -gt 0 ]]
	then 
		THREEDDNA_VISUALIZE_OPT="${THREEDDNA_VISUALIZE_OPT} -i"
	fi	
}    	

if [[ -z ${SALSA_PATH} || ! -f ${SALSA_PATH}/run_pipeline.py ]]
then
	(>&2 echo "Variable SALSA_PATH must be set to a proper salsa2 installation directory!!")
    exit 1
fi

if [[ -z ${THREEDDNA_PATH} || ! -f ${THREEDDNA_PATH}/run-asm-pipeline.sh ]]
then 
	(>&2 echo "Variable THREEDDNA_PATH must be set to a proper 3d-dna installation directory!!")
    exit 1
fi
 
# Type: 0 - Arima Mapping Pipeline (For QC) + SALSA SCAFFOLDING 
# Type: 1 - Phase Genomics Mapping Pipeline (For QC)
# Type: 2 - Aiden Lab Juicer/3d-dna Pipeline (For QC) - full scaffolding pipeline
# Type: 3 - Aiden Lab fast Juicer + visualize any input assembly 
# Type: 4 - higlass visualization

myTypes=("01_HICsalsaPrepareInput, 02_HICsalsaBwa, 03_HICsalsaFilter, 04_HICsalsaMerge, 05_HICsalsaMarkduplicates, 06_HICsalsaSalsa, 07_HICsalsaStatistics", 
"01_HICphasePrepareInput, 02_HICphaseBwa, 03_HICphaseFilter, 04_HICphaseMatlock", 
"01_HIC3dnaPrepareInput, 02_HIC3dnaJuicer, 03_HIC3dnaAssemblyPipeline",
"01_HIC3dnaPrepareInput, 02_HIC3dnaJuicer, 03_HIC3dnaVisualize",
"01_HIChiglassPrepare, 02_HiChiglassBwa, 03_HiChiglassFilter, 04_HiChiglassMatrix")
if [[ ${SC_HIC_TYPE} -eq 0 ]]
then 
    ### 01_HICsalsaPrepareInput
    if [[ ${currentStep} -eq 1 ]]
    then
        ### clean up plans 
        for x in $(ls hic_01_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        if [[ ! -f "${SC_HIC_REF}" ]]
        then
        	(>&2 echo "ERROR - set SC_HIC_REF to reference fasta file")
        	exit 1
   		fi
   		
   		if [[ ! -d "${SC_HIC_READS}" ]]
        then
        	(>&2 echo "ERROR - set SC_HIC_READS to HiC read directory")
        	exit 1
   		fi
   		
   		echo "if [[ -d ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID} ]]; then mv ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID} ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}_$(date '+%Y-%m-%d_%H-%M-%S'); fi && mkdir ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}" > hic_01_HICsalsaPrepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "mkdir -p ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/reads" >> hic_01_HICsalsaPrepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "mkdir -p ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/bams" >> hic_01_HICsalsaPrepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "mkdir -p ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/ref" >> hic_01_HICsalsaPrepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "mkdir -p ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/config" >> hic_01_HICsalsaPrepareInput_single_${CONT_DB}.${slurmID}.plan
		
		if [[ -f ${SC_HIC_REF_EXCLUDELIST} ]]
   		then
   			echo "${SEQKIT_PATH} grep -v -f ${SC_HIC_REF_EXCLUDELIST} ${SC_HIC_REF} | sed -e \"s/:/-/g\" ${SC_HIC_REF} > ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/ref/$(basename ${SC_HIC_REF})" 
   			echo "${SEQKIT_PATH} grep -f ${SC_HIC_REF_EXCLUDELIST} ${SC_HIC_REF} | sed -e \"s/:/-/g\" ${SC_HIC_REF} > ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/ref/exclude.fasta" 
   		else   			
   			echo "sed -e \"s/:/-/g\" ${SC_HIC_REF} > ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/ref/$(basename ${SC_HIC_REF})"
   		fi >> hic_01_HICsalsaPrepareInput_single_${CONT_DB}.${slurmID}.plan
   		
   		if [[ -v ${SC_HIC_REF_FIXBIONANOGAPS} && ${SC_HIC_REF_FIXBIONANOGAPS} -gt 0 ]]
   		then
   			tmpName=${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/ref/$(basename ${SC_HIC_REF%.fasta})
   			tmpName=${tmpName%.fa}_tmp.fasta
   			echo "mv ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/ref/$(basename ${SC_HIC_REF}) ${tmpName}"
   			echo "${SEQKIT_PATH} replace -p \"^[Nn]+|[nN]+$\" -r \"\" -s ${tmpName} | ${SEQKIT_PATH} replace -p \"[Nn]+\" -r \"\$(printf 'N%.0s' {1..111})\" -s > ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/ref/$(basename ${SC_HIC_REF})"   			
   		fi >> hic_01_HICsalsaPrepareInput_single_${CONT_DB}.${slurmID}.plan
   		
		echo "samtools faidx ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/ref/$(basename ${SC_HIC_REF})" >> hic_01_HICsalsaPrepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "bwa index ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/ref/$(basename ${SC_HIC_REF})" >> hic_01_HICsalsaPrepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "cp ${configFile} ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/config/$(basename ${configFile%.sh})_$(date '+%Y-%m-%d_%H-%M-%S').sh" >> hic_01_HICsalsaPrepareInput_single_${CONT_DB}.${slurmID}.plan
		
		echo "samtools $(${CONDA_HIC_ENV} && samtools 2>&1 | grep Version | awk '{print $2}' && conda deactivate)" > hic_01_HICsalsaPrepareInput_single_${CONT_DB}.${slurmID}.version
		echo "bwa $(${CONDA_HIC_ENV} && bwa 2>&1 | grep Version | awk '{print $2}' && conda deactivate)" >> hic_01_HICsalsaPrepareInput_single_${CONT_DB}.${slurmID}.version
	### 02_HICsalsaBwa 
    elif [[ ${currentStep} -eq 2 ]]
    then
        ### clean up plans 
        for x in $(ls hic_02_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
		
		if [[ ! -d "${SC_HIC_READS}" ]]
        then
        	(>&2 echo "ERROR - set SC_HIC_READS to directory that contain the PROJECT_ID*.fastq.qz read files")
        	exit 1
   		fi   		   				
        
        if [[ ! -d "${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/reads" ]]
        then
        	(>&2 echo "ERROR - cannot access directory ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/reads!")
        	exit 1
   		fi
   		
   		ref=${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/ref/$(basename ${SC_HIC_REF})
   		
   		if [[ ! -f "${ref}" ]]
        then
        (>&2 echo "ERROR - cannot reference fasta file ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/ref/$(basename ${SC_HIC_REF})!")
        	exit 1
   		fi
   		
   		### link HiC reads into current reads sub directory
   		for x in ${SC_HIC_READS}/${PROJECT_ID}_*_*_R[12].fastq.gz
		do
			if [[ -f ${x} ]]
			then	
				ln -s -r -f ${x} ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/reads 
			fi
		done
   		   				
   		numR1Files=0
		for x in ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/reads/${PROJECT_ID}_*_*_R1.fastq.gz
		do
			if [[ -f ${x} ]]
			then	
				numR1Files=$((${numR1Files}+1))	
			fi
		done
		
		if [[ ${numR1Files} -eq 0 ]]
        then
        	(>&2 echo "ERROR - cannot read HiC R1 files with following pattern: ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/reads/${PROJECT_ID}_*_*_R1.fastq.gz")
        	exit 1
   		fi
   		
   		numR2Files=0
		for x in ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/reads/${PROJECT_ID}_*_*_R2.fastq.gz
		do
			if [[ -f ${x} ]]
			then	
				numR2Files=$((${numR2Files}+1))	
			fi
		done
		
		if [[ ${numR2Files} -eq 0 ]]
        then
        	(>&2 echo "ERROR - cannot read HiC R2 files with following pattern: ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/reads/${PROJECT_ID}_*_*_R2.fastq.gz")
        	exit 1
   		fi
   		
   		if [[ ${numR1Files} -ne ${numR2Files} ]]
        then
        	(>&2 echo "ERROR - HiC R1 files ${numR1Files} does not match R2 files ${numR2Files}")
        	exit 1
   		fi
   		
   		setbwaOptions
   		
		for r1 in ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/reads/${PROJECT_ID}_*_*_R1.fastq.gz
		do
			id=$(dirname ${r1})
			f1=$(basename ${r1})
			f2=$(echo "${f1}" | sed -e "s:_R1.fastq.gz:_R2.fastq.gz:")
			o="${f1%_R1.fastq.gz}"											
			
			echo "bwa mem${CONTIG_BWA_OPT} -R \"@RG\tID:${o}\tSM:${PROJECT_ID}_HIC\tLB:${PROJECT_ID}_HIC\tPL:ILLUMINA\tPU:none\" ${ref} ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/reads/${f1} | samtools view -Sb - > ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/bams/${o}_bwa_1.bam"
			echo "bwa mem${CONTIG_BWA_OPT} -R \"@RG\tID:${o}\tSM:${PROJECT_ID}_HIC\tLB:${PROJECT_ID}_HIC\tPL:ILLUMINA\tPU:none\" ${ref} ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/reads/${f2} | samtools view -Sb - > ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/bams/${o}_bwa_2.bam" 				 
		done > hic_02_HICsalsaBwa_block_${CONT_DB}.${slurmID}.plan
		
   		echo "bwa $(${CONDA_HIC_ENV} && bwa 2>&1 | grep Version | awk '{print $2}' && conda deactivate)" > hic_02_HICsalsaBwa_block_${CONT_DB}.${slurmID}.version
   		echo "samtools $(${CONDA_HIC_ENV} && samtools 2>&1 | grep Version | awk '{print $2}' && conda deactivate)" >> hic_02_HICsalsaBwa_block_${CONT_DB}.${slurmID}.version
	### 03_HICsalsaFilter
    elif [[ ${currentStep} -eq 3 ]]
    then
        ### clean up plans 
        for x in $(ls hic_03_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done

		if [[ ! -d "${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/bams" ]]
        then
        	(>&2 echo "ERROR - cannot access directory ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/bams!")
        	exit 1
   		fi
   		   		   				
		for b1 in ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/bams/*_bwa_1.bam
		do
			d=$(dirname ${b1})
			b2="${b1%_1.bam}_2.bam"
			f1=$(basename ${b1%_bwa_1.bam})_bwaFilt_1.bam
			f2="${f1%_1.bam}_2.bam"			
			
			echo "samtools view -h ${b1} | perl ${MARVEL_PATH}/scripts/filter_five_end.pl | samtools view -Sb - > ${d}/${f1}"
			echo "samtools view -h ${b2} | perl ${MARVEL_PATH}/scripts/filter_five_end.pl | samtools view -Sb - > ${d}/${f2}" 				 
		done > hic_03_HICsalsaFilter_block_${CONT_DB}.${slurmID}.plan
					
		echo "samtools $(${CONDA_HIC_ENV} && samtools 2>&1 | grep Version | awk '{print $2}' && conda deactivate)" > hic_03_HICsalsaFilter_block_${CONT_DB}.${slurmID}.version	   	
	### 04_HICsalsaMerge
    elif [[ ${currentStep} -eq 4 ]]
    then
        ### clean up plans 
        for x in $(ls hic_04_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        if [[ ! -d "${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/bams" ]]
        then
        	(>&2 echo "ERROR - cannot access directory ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/bams!")
        	exit 1
   		fi
   		
   		ref=${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/ref/$(basename ${SC_HIC_REF})
   		
   		if [[ ! -f ${ref}.fai ]]
   		then  
   		 	(>&2 echo "ERROR - cannot access reference fasta index ${ref}.fai!")
        	exit 1
		fi
		 
		if [[ -z ${SC_HIC_MINMAPQV} ]]
		then
			SC_HIC_MINMAPQV=10	
		fi
		 
		for b1 in ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/bams/*_bwaFilt_1.bam
		do
			b2="${b1%_1.bam}_2.bam"
			o="${b1%_1.bam}.bam"			
			
			echo "perl ${MARVEL_PATH}/scripts/two_read_bam_combiner.pl ${b1} ${b2} $(which samtools) ${SC_HIC_MINMAPQV} | samtools view -bS -t ${ref}.fai - | samtools sort -o ${o} -"			 				 
			done > hic_04_HICsalsaMerge_single_${CONT_DB}.${slurmID}.plan
		   
		echo "samtools $(${CONDA_HIC_ENV} && samtools 2>&1 | grep Version | awk '{print $2}' && conda deactivate)" > hic_04_HICsalsaMerge_single_${CONT_DB}.${slurmID}.version		
   	### 05_HICsalsaMarkduplicates
    elif [[ ${currentStep} -eq 5 ]]
    then
        ### clean up plans 
        for x in $(ls hic_05_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
		if [[ ! -d "${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/bams" ]]
        then
        	(>&2 echo "ERROR - cannot access directory ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/bams!")
        	exit 1
   		fi
   		
   		setPicardOptions
   		   		
   		## if multiple bam files are available (e.g. different Lanes) then merge files prior to markduplicates
   		files=$(ls ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/bams/*_bwaFilt.bam)
   		   		   		
   		echo "picard MarkDuplicates $(${CONDA_HIC_ENV} && picard MarkDuplicates --version && conda deactivate)" > hic_05_HICsalsaMarkduplicates_single_${CONT_DB}.${slurmID}.version
   		
   		if [[ $(echo $files | wc -w) -eq 1 ]]
   		then
   			ob="${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/bams/${PROJECT_ID}_finalHiC.bam"
			m="${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/bams/${PROJECT_ID}_finalHiC.metrics"
		   echo "picard ${CONTIG_PICARD_OPT} MarkDuplicates I=${files} O=${ob} M=${m} && samtools index -@ ${SC_HIC_SAMTOOLS_THREADS} ${ob} && ln -s -f -r ${ob} ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID} && ln -s -f -r ${ob}.bai ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}"
   		elif [[ $(echo $files | wc -w) -gt 1 ]]
   		then
   			mrg=${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/bams/${PROJECT_ID}_mergedHiC.bam
   			o=${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/bams/${PROJECT_ID}_finalHiC.bam
   			m=${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/bams/${PROJECT_ID}_finalHiC.metrics
   			i=$(echo -e ${files} | sed -e "s:${SC_HIC_OUTDIR}/:I=${SC_HIC_OUTDIR}/:g")
   			echo "picard ${CONTIG_PICARD_OPT} MergeSamFiles ${i} OUTPUT=${mrg} USE_THREADING=TRUE ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT && picard ${CONTIG_PICARD_OPT} MarkDuplicates I=${mrg} O=${o} M=${m} && samtools index -@ ${SC_HIC_SAMTOOLS_THREADS} ${o} && ln -s -f -r ${o} ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID} && ln -s -f -r ${ob}.bai ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}"
   			echo "picard MergeSamFiles $(${CONDA_HIC_ENV} && picard MarkDuplicates --version && conda deactivate)" >> hic_05_HICsalsaMarkduplicates_single_${CONT_DB}.${slurmID}.version	
   		else
   	 		(>&2 echo "ERROR - cannot find file with following pattern: ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/bams/*_bwaFilt.bam!")
        	exit 1
   	 	fi > hic_05_HICsalsaMarkduplicates_single_${CONT_DB}.${slurmID}.plan 
   	### 06_HICsalsaSalsa
    elif [[ ${currentStep} -eq 6 ]]
    then
        ### clean up plans 
        for x in $(ls hic_06_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
      	done
       
       	if [[ ! -f "${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/bams/${PROJECT_ID}_finalHiC.bam" ]]
       	then
    		(>&2 echo "ERROR - cannot access final duplicate marked bam file ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/bams/${PROJECT_ID}_finalHiC.bam!")
        	exit 1
		fi
		
		ref="${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/ref/$(basename ${SC_HIC_REF})"
		
		if [[ ! -f ${ref} ]]
       	then
    		(>&2 echo "ERROR - cannot access reference fasta file: \"${ref}\"!")
        	exit 1
		fi
		
		if [[ ! -f ${ref}.fai ]]
       	then
    		(>&2 echo "ERROR - cannot access reference fasta index file: \"${ref}.fai\"!")
        	exit 1
		fi
		
		if [[ -z ${SC_HIC_ENZYME_SEQ} ]]
       	then 
       		(>&2 echo "ERROR - Enzyme is required, set variable SC_HIC_ENZYME_SEQ!")
        	exit 1	
       	fi
       	
       	SALSA_OPT="-s ${gsize}"
       	if [[ -v ${SC_HIC_REF_DOBREAKS} && ${SC_HIC_REF_DOBREAKS} -gt 0 ]]
       	then
       		SALSA_OPT="${SALSA_OPT} -m yes"        		
       	fi
       	
       	if [[ -v ${SC_HIC_REF_NUMITER} && ${SC_HIC_REF_NUMITER} -gt 0 ]]
       	then
       		SALSA_OPT="${SALSA_OPT} -i ${SC_HIC_REF_NUMITER}"     		
       	fi
       	
		if [[ -v ${SC_HIC_MINCONTIG} && ${SC_HIC_MINCONTIG} -gt 0 ]]
       	then
       		SALSA_OPT="${SALSA_OPT} -c ${SC_HIC_MINCONTIG}"     		
       	fi
       	
       	if [[ -v ${SC_HIC_KEEPINTERMEDIATES} && ${SC_HIC_KEEPINTERMEDIATES} -gt 0 ]]
       	then
       		SALSA_OPT="${SALSA_OPT} -p yes"     		
       	fi
		
		echo "bedtools bamtobed -i ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/bams/${PROJECT_ID}_finalHiC.bam > ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/bams/${PROJECT_ID}_finalHiC.bed" > hic_06_HICsalsaSalsa_single_${CONT_DB}.${slurmID}.plan
		echo "sort -k 4 ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/bams/${PROJECT_ID}_finalHiC.bed > ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/bams/${PROJECT_ID}_finalHiC_sortByName.bed" >> hic_06_HICsalsaSalsa_single_${CONT_DB}.${slurmID}.plan
       	bed=${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/bams/${PROJECT_ID}_finalHiC_sortByName.bed
       	echo "export PATH=${SALSA_PATH}:\$PATH && run_pipeline.py -a ${ref} -l ${ref}.fai -b ${bed} -e ${SC_HIC_ENZYME_SEQ} -o ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/out ${SALSA_OPT}" >> hic_06_HICsalsaSalsa_single_${CONT_DB}.${slurmID}.plan
       	
       	echo "SALSA $(git --git-dir=${SALSA_PATH}/.git rev-parse --short HEAD)" > hic_06_HICsalsaSalsa_single_${CONT_DB}.${slurmID}.version
       	echo "${CONDA_HIC_ENV} && bedtools --version && conda deactivate" >> hic_06_HICsalsaSalsa_single_${CONT_DB}.${slurmID}.version
    ### 07_HICsalsaStatistics
    elif [[ ${currentStep} -eq 7 ]]
    then
        ### clean up plans 
        for x in $(ls hic_07_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
		if [[ ! -f "${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/bams/${PROJECT_ID}_finalHiC.bam" ]]
        then
    	(>&2 echo "ERROR - cannot access final duplicate marked bam file ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/bams/${PROJECT_ID}_finalHiC.bam!")
        	exit 1
   		fi
   		   		  
   		echo "perl ${MARVEL_PATH}/scripts/get_stats.pl ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/bams/${PROJECT_ID}_finalHiC.bam > ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/${PROJECT_ID}_finalHiC.stats" > hic_07_HICsalsaStatistics_single_${CONT_DB}.${slurmID}.plan
   		
   		if [[ -n ${SC_HIC_FULLSTATS} && ${SC_HIC_FULLSTATS} -gt 0 ]]
   		then
			### run slurm stats - on the master node !!! Because sacct is not available on compute nodes
	    	if [[ $(hostname) == "falcon1" || $(hostname) == "falcon2" ]]
	        then 
	        	bash ${SUBMIT_SCRIPTS_PATH}/slurmStats.sh ${configFile}
	    	else
	        	cwd=$(pwd)
	        	ssh falcon "cd ${cwd} && bash ${SUBMIT_SCRIPTS_PATH}/slurmStats.sh ${configFile}"
	    	fi
		fi
    	### create assemblyStats plan 
    	echo "${SUBMIT_SCRIPTS_PATH}/assemblyStats.sh ${configFile} 14" >> hic_07_HICsalsaStatistics_single_${CONT_DB}.${slurmID}.plan
    	git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD >> hic_07_HICsalsaStatistics_single_${CONT_DB}.${slurmID}.version   		  	
    else
        (>&2 echo "step ${currentStep} in SC_HIC_TYPE ${SC_HIC_TYPE} not supported")
        (>&2 echo "valid steps are: ${myTypes[${SC_HIC_TYPE}]}")
        exit 1            
    fi 
elif [[ ${SC_HIC_TYPE} -eq 1 ]] ### 01_HICphasePrepareInput, 02_HICphaseBwa, 03_HICphaseFilter, 04_HICphaseMatlock
then     
 	(>&2 echo "Phase qc not implemented yet")
    exit 1
elif [[ ${SC_HIC_TYPE} -eq 2 ]] 
then 
	### 01_HIC3dnaPrepareInput
	if [[ ${currentStep} -eq 1 ]]
    then
        ### clean up plans 
        for x in $(ls hic_01_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        if [[ -z ${JUICER_PATH} || ! -f ${JUICER_PATH}/SLURM/scripts/juicer.sh || -z ${JUICER_TOOLS_PATH} || ! -f ${JUICER_TOOLS_PATH} ]]
        then 
    		(>&2 echo "[ERROR] Set variable JUICER_PATH to juicer install directory and JUICER_TOOLS_PATH to a valid jar file (e.g. juicer_tools.1.9.8_jcuda.0.8.jar)")
    		exit 1
    	fi
    	
    	
    	if [[ ! -f "${SC_HIC_REF}" ]]
        then
        	(>&2 echo "ERROR - set SC_HIC_REF to reference fasta file")
        	exit 1
   		fi
   		
   		if [[ ! -d "${SC_HIC_READS}" ]]
        then
        	(>&2 echo "ERROR - set SC_HIC_READS to HiC read directory")
        	exit 1
   		fi
   		
		numR1Files=0
		for x in ${SC_HIC_READS}/${PROJECT_ID}_*_*_R1.fastq.gz
		do
			if [[ -f ${x} ]]
			then	
				numR1Files=$((${numR1Files}+1))	
			fi
		done
		
		if [[ ${numR1Files} -eq 0 ]]
        then
        	(>&2 echo "ERROR - cannot read HiC R1 files with following pattern: ${SC_HIC_READS}/${PROJECT_ID}_*_*_R1.fastq.gz")
        	exit 1
   		fi
		
		numR2Files=0
		for x in ${SC_HIC_READS}/${PROJECT_ID}_*_*_R2.fastq.gz
		do
			if [[ -f ${x} ]]
			then	
				numR2Files=$((${numR2Files}+1))	
			fi
		done
		
		if [[ ${numR2Files} -eq 0 ]]
        then
        	(>&2 echo "ERROR - cannot read HiC R2 files with following pattern: ${SC_HIC_READS}/${PROJECT_ID}_*_*_R2.fastq.gz")
        	exit 1
   		fi
   		
   		if [[ ${numR1Files} -ne ${numR2Files} ]]
        then
        	(>&2 echo "ERROR - HiC R1 files ${numR1Files} does not match R2 files ${numR2Files}")
        	exit 1
   		fi   		
   		
   		if [[ -z ${SC_HIC_ENZYME_NAME} ]]
   		then
        	(>&2 echo "ERROR - set variable SC_HIC_ENZYME_NAME!")
        	exit 1   			
   		fi
   		
   		echo "if [[ -d ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID} ]]; then mv ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID} ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}_$(date '+%Y-%m-%d_%H-%M-%S'); fi && mkdir ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}" > hic_01_HIC3dnaPrepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "cp -r ${JUICER_PATH}/SLURM/scripts ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}" >> hic_01_HIC3dnaPrepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "ln -s -f -r ${JUICER_TOOLS_PATH} ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/scripts/juicer_tools.jar" >> hic_01_HIC3dnaPrepareInput_single_${CONT_DB}.${slurmID}.plan
		# create reference subdir + link and do indexing 
		echo "mkdir -p ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/references" >> hic_01_HIC3dnaPrepareInput_single_${CONT_DB}.${slurmID}.plan
		
		if [[ -f ${SC_HIC_REF_EXCLUDELIST} ]]
   		then
   			echo "${SEQKIT_PATH} grep -v -f ${SC_HIC_REF_EXCLUDELIST} ${SC_HIC_REF} > ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/references/${PROJECT_ID}.fasta" >>hic_01_HIC3dnaPrepareInput_single_${CONT_DB}.${slurmID}.plan
   			echo "${SEQKIT_PATH} grep -f ${SC_HIC_REF_EXCLUDELIST} ${SC_HIC_REF} > ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/references/exclude.fasta" >> hic_01_HIC3dnaPrepareInput_single_${CONT_DB}.${slurmID}.plan
   		else   			
   			echo "ln -s -r ${SC_HIC_REF} ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/references/${PROJECT_ID}.fasta" >> hic_01_HIC3dnaPrepareInput_single_${CONT_DB}.${slurmID}.plan
   		fi
		
		echo "samtools faidx ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/references/${PROJECT_ID}.fasta" >> hic_01_HIC3dnaPrepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "awk '{print \$1\" \"\$2}' ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/references/${PROJECT_ID}.fasta.fai > ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/references/${PROJECT_ID}.sizes" >> hic_01_HIC3dnaPrepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "bwa index ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/references/${PROJECT_ID}.fasta" >> hic_01_HIC3dnaPrepareInput_single_${CONT_DB}.${slurmID}.plan
		# create resctriction site subdir + generate sites
		echo "mkdir -p ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/restriction_sites" >> hic_01_HIC3dnaPrepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "python2 ${MARVEL_PATH}/scripts/generate_site_positions.py ${SC_HIC_ENZYME_NAME} ${PROJECT_ID} ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/references/${PROJECT_ID}.fasta ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/restriction_sites" >> hic_01_HIC3dnaPrepareInput_single_${CONT_DB}.${slurmID}.plan
		# create fastq subdir + link zipped HIC reads  
		echo "mkdir -p ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/fastq" >> hic_01_HIC3dnaPrepareInput_single_${CONT_DB}.${slurmID}.plan
   		for x in ${SC_HIC_READS}/${PROJECT_ID}_*_*_R[12].fastq.gz
		do
			if [[ -f ${x} ]]
			then	
				echo "ln -s -r -f ${x} ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/fastq" >> hic_01_HIC3dnaPrepareInput_single_${CONT_DB}.${slurmID}.plan 
			fi
		done		
		
		# create config subdir + copy current config file with time stamp
		echo "mkdir -p ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/config" >> hic_01_HIC3dnaPrepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "cp ${configFile} ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/config/$(basename ${configFile%.sh})_$(date '+%Y-%m-%d_%H-%M-%S').sh" >> hic_01_HIC3dnaPrepareInput_single_${CONT_DB}.${slurmID}.plan
		
		echo "samtools $(${CONDA_HIC_ENV} && samtools 2>&1 | grep Version | awk '{print $2}' && conda deactivate)" > hic_01_HIC3dnaPrepareInput_single_${CONT_DB}.${slurmID}.version
		echo "bwa $(${CONDA_HIC_ENV} && bwa 2>&1 | grep Version | awk '{print $2}' && conda deactivate)" >> hic_01_HIC3dnaPrepareInput_single_${CONT_DB}.${slurmID}.version
	### 02_HIC3dnaJuicer
	elif [[ ${currentStep} -eq 2 ]]
    then
        ### clean up plans 
        for x in $(ls hic_02_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        setJuicerOptions
        echo "${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/scripts/juicer.sh ${SC_HIC_JUICER_OPT} -d $(pwd)/${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID} -D $(pwd)/${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID} -g ${PROJECT_ID} -s ${SC_HIC_ENZYME_NAME} -z ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/references/${PROJECT_ID}.fasta -y ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/restriction_sites/${PROJECT_ID}_${SC_HIC_ENZYME_NAME}.txt -p ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/references/${PROJECT_ID}.sizes" > hic_02_HIC3dnaJuicer_single_${CONT_DB}.${slurmID}.plan    	        
    ### 03_HIC3dnaAssemblyPipeline
	elif [[ ${currentStep} -eq 3 ]]
    then
        ### clean up plans 
        for x in $(ls hic_03_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        setThreeDDNAOptions
            	        
    	echo "cd ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID} && ${THREEDDNA_PATH}/run-asm-pipeline.sh${THREEDDNA_OPT} references/${PROJECT_ID}.fasta aligned/merged_nodups.txt" > hic_03_HIC3dnaAssemblyPipeline_single_${CONT_DB}.${slurmID}.plan
        
        echo "3d-dna $(git --git-dir=${THREEDDNA_PATH}/.git rev-parse --short HEAD)" > hic_03_HIC3dnaAssemblyPipeline_single_${CONT_DB}.${slurmID}.version
  	else
    	(>&2 echo "step ${currentStep} in SC_HIC_TYPE ${SC_HIC_TYPE} not supported")
    	(>&2 echo "valid steps are: ${myTypes[${SC_HIC_TYPE}]}")
    	exit 1
	fi
elif [[ ${SC_HIC_TYPE} -eq 3 ]] 
then 
	### 01_HIC3dnaPrepareInput
	if [[ ${currentStep} -eq 1 ]]
    then
        ### clean up plans 
        for x in $(ls hic_01_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        if [[ -z ${JUICER_PATH} || ! -f ${JUICER_PATH}/SLURM/scripts/juicer.sh || -z ${JUICER_TOOLS_PATH} || ! -f ${JUICER_TOOLS_PATH} ]]
        then 
    		(>&2 echo "[ERROR] Set variable JUICER_PATH to juicer install directory and JUICER_TOOLS_PATH to a valid jar file (e.g. juicer_tools.1.9.8_jcuda.0.8.jar)")
    		exit 1
    	fi
    	
    	
    	if [[ ! -f "${SC_HIC_REF}" ]]
        then
        	(>&2 echo "ERROR - set SC_HIC_REF to reference fasta file")
        	exit 1
   		fi
   		
   		if [[ ! -d "${SC_HIC_READS}" ]]
        then
        	(>&2 echo "ERROR - set SC_HIC_READS to HiC read directory")
        	exit 1
   		fi
   		
		numR1Files=0
		for x in ${SC_HIC_READS}/${PROJECT_ID}_*_*_R1.fastq.gz
		do
			if [[ -f ${x} ]]
			then	
				numR1Files=$((${numR1Files}+1))	
			fi
		done
		
		if [[ ${numR1Files} -eq 0 ]]
        then
        	(>&2 echo "ERROR - cannot read HiC R1 files with following pattern: ${SC_HIC_READS}/${PROJECT_ID}_*_*_R1.fastq.gz")
        	exit 1
   		fi
		
		numR2Files=0
		for x in ${SC_HIC_READS}/${PROJECT_ID}_*_*_R2.fastq.gz
		do
			if [[ -f ${x} ]]
			then	
				numR2Files=$((${numR2Files}+1))	
			fi
		done
		
		if [[ ${numR2Files} -eq 0 ]]
        then
        	(>&2 echo "ERROR - cannot read HiC R2 files with following pattern: ${SC_HIC_READS}/${PROJECT_ID}_*_*_R2.fastq.gz")
        	exit 1
   		fi
   		
   		if [[ ${numR1Files} -ne ${numR2Files} ]]
        then
        	(>&2 echo "ERROR - HiC R1 files ${numR1Files} does not match R2 files ${numR2Files}")
        	exit 1
   		fi   		
   		
   		if [[ -z ${SC_HIC_ENZYME_NAME} ]]
   		then
        	(>&2 echo "ERROR - set variable SC_HIC_ENZYME_NAME!")
        	exit 1   			
   		fi
   		
   		echo "if [[ -d ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID} ]]; then mv ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID} ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}_$(date '+%Y-%m-%d_%H-%M-%S'); fi && mkdir ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}" > hic_01_HIC3dnaPrepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "cp -r ${JUICER_PATH}/SLURM/scripts ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}" >> hic_01_HIC3dnaPrepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "ln -s -f -r ${JUICER_TOOLS_PATH} ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/scripts/juicer_tools.jar" >> hic_01_HIC3dnaPrepareInput_single_${CONT_DB}.${slurmID}.plan
		# create reference subdir + link and do indexing 
		echo "mkdir -p ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/references" >> hic_01_HIC3dnaPrepareInput_single_${CONT_DB}.${slurmID}.plan
		
		if [[ -f ${SC_HIC_REF_EXCLUDELIST} ]]
   		then
   			echo "${SEQKIT_PATH} grep -v -f ${SC_HIC_REF_EXCLUDELIST} ${SC_HIC_REF} > ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/references/${PROJECT_ID}.fasta" >>hic_01_HIC3dnaPrepareInput_single_${CONT_DB}.${slurmID}.plan
   			echo "${SEQKIT_PATH} grep -f ${SC_HIC_REF_EXCLUDELIST} ${SC_HIC_REF} > ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/references/exclude.fasta" >> hic_01_HIC3dnaPrepareInput_single_${CONT_DB}.${slurmID}.plan
   		else   			
   			echo "ln -s -r ${SC_HIC_REF} ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/references/${PROJECT_ID}.fasta" >> hic_01_HIC3dnaPrepareInput_single_${CONT_DB}.${slurmID}.plan
   		fi		
		
		echo "samtools faidx ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/references/${PROJECT_ID}.fasta" >> hic_01_HIC3dnaPrepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "awk '{print \$1\" \"\$2}' ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/references/${PROJECT_ID}.fasta.fai > ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/references/${PROJECT_ID}.sizes" >> hic_01_HIC3dnaPrepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "bwa index ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/references/${PROJECT_ID}.fasta" >> hic_01_HIC3dnaPrepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "awk -f ${THREEDDNA_PATH}/utils/generate-assembly-file-from-fasta.awk ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/references/${PROJECT_ID}.fasta > ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/references/${PROJECT_ID}.assembly" >> hic_01_HIC3dnaPrepareInput_single_${CONT_DB}.${slurmID}.plan
		# create resctriction site subdir + generate sites
		echo "mkdir -p ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/restriction_sites" >> hic_01_HIC3dnaPrepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "python ${MARVEL_PATH}/scripts/generate_site_positions.py ${SC_HIC_ENZYME_NAME} ${PROJECT_ID} ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/references/${PROJECT_ID}.fasta ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/restriction_sites" >> hic_01_HIC3dnaPrepareInput_single_${CONT_DB}.${slurmID}.plan
		# create fastq subdir + link zipped HIC reads  
		echo "mkdir -p ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/fastq" >> hic_01_HIC3dnaPrepareInput_single_${CONT_DB}.${slurmID}.plan
   		for x in ${SC_HIC_READS}/${PROJECT_ID}_*_*_R[12].fastq.gz
		do
			if [[ -f ${x} ]]
			then	
				echo "ln -s -r -f ${x} ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/fastq" >> hic_01_HIC3dnaPrepareInput_single_${CONT_DB}.${slurmID}.plan 
			fi
		done
		echo "mkdir -p ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/visualize" >> hic_01_HIC3dnaPrepareInput_single_${CONT_DB}.${slurmID}.plan		
		
		# create config subdir + copy current config file with time stamp
		echo "mkdir -p ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/config" >> hic_01_HIC3dnaPrepareInput_single_${CONT_DB}.${slurmID}.plan
		echo "cp ${configFile} ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/config/$(basename ${configFile%.sh})_$(date '+%Y-%m-%d_%H-%M-%S').sh" >> hic_01_HIC3dnaPrepareInput_single_${CONT_DB}.${slurmID}.plan
		
		echo "samtools $(${CONDA_HIC_ENV} && samtools 2>&1 | grep Version | awk '{print $2}' && conda deactivate)" > hic_01_HIC3dnaPrepareInput_single_${CONT_DB}.${slurmID}.version
		echo "bwa $(${CONDA_HIC_ENV} && bwa 2>&1 | grep Version | awk '{print $2}' && conda deactivate)" >> hic_01_HIC3dnaPrepareInput_single_${CONT_DB}.${slurmID}.version
		echo "3d-dna $(git --git-dir=${THREEDDNA_PATH}/.git rev-parse --short HEAD)" >> hic_01_HIC3dnaPrepareInput_single_${CONT_DB}.${slurmID}.version
	### 02_HIC3dnaJuicer
	elif [[ ${currentStep} -eq 2 ]]
    then
        ### clean up plans 
        for x in $(ls hic_02_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        setJuicerOptions 1
        
        if [[ -d ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/aligned ]]
        then
        	mv ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/aligned ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/aligned_old_$(date '+%Y-%m-%d_%H-%M-%S')
    	fi    	        
            	        
        echo "${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/scripts/juicer.sh ${SC_HIC_JUICER_OPT} -S early -d $(pwd)/${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID} -D $(pwd)/${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID} -g ${PROJECT_ID} -s ${SC_HIC_ENZYME_NAME} -z ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/references/${PROJECT_ID}.fasta -y ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/restriction_sites/${PROJECT_ID}_${SC_HIC_ENZYME_NAME}.txt -p ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/references/${PROJECT_ID}.sizes" > hic_02_HIC3dnaJuicer_single_${CONT_DB}.${slurmID}.plan
    ### 03_HIC3dnaVisualize
	elif [[ ${currentStep} -eq 3 ]]
    then
        ### clean up plans 
        for x in $(ls hic_03_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
            
        set3DDNAVisualizeOptions  
        echo "mkdir -p ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/visualize" > hic_03_HIC3dnaVisualizePipeline_single_${CONT_DB}.${slurmID}.plan    	        
        if [[ ! -f ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/references/${PROJECT_ID}.assembly ]]
        then
        	echo "awk -f ${THREEDDNA_PATH}/utils/generate-assembly-file-from-fasta.awk ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/references/${PROJECT_ID}.fasta > ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/references/${PROJECT_ID}.assembly" >> hic_03_HIC3dnaVisualizePipeline_single_${CONT_DB}.${slurmID}.plan
    	fi
    	echo "cd ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/visualize && ${THREEDDNA_PATH}/visualize/run-assembly-visualizer.sh ${THREEDDNA_VISUALIZE_OPT} ../references/${PROJECT_ID}.assembly ../aligned/merged_nodups.txt" >> hic_03_HIC3dnaVisualizePipeline_single_${CONT_DB}.${slurmID}.plan
        
        echo "3d-dna $(git --git-dir=${THREEDDNA_PATH}/.git rev-parse --short HEAD)" > hic_03_HIC3dnaVisualizePipeline_single_${CONT_DB}.${slurmID}.version
  	else
    	(>&2 echo "step ${currentStep} in SC_HIC_TYPE ${SC_HIC_TYPE} not supported")
    	(>&2 echo "valid steps are: ${myTypes[${SC_HIC_TYPE}]}")
    	exit 1
	fi	
### 01_HIChiglassPrepare, 02_HiChiglassBwa, 03_HiChiglassFilter, 04_HiChiglassMatrix	
elif [[ ${SC_HIC_TYPE} -eq 4 ]] 
then 
	### 01_HIChiglassPrepare, 02_HiChiglassBwa, 03_HiChiglassFilter, 04_HiChiglassMatrix
	if [[ ${currentStep} -eq 1 ]]
    then
        ### clean up plans 
        for x in $(ls hic_01_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        if [[ ! -f "${SC_HIC_REF}" ]]
        then
        	(>&2 echo "ERROR - set SC_HIC_REF to reference fasta file")
        	exit 1
   		fi
   		
   		if [[ ! -d "${SC_HIC_READS}" ]]
        then
        	(>&2 echo "ERROR - set SC_HIC_READS to HiC read directory")
        	exit 1
   		fi
   		
   		echo "if [[ -d ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID} ]]; then mv ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID} ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}_$(date '+%Y-%m-%d_%H-%M-%S'); fi && mkdir ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}" > hic_01_HIChiglassPrepare_single_${CONT_DB}.${slurmID}.plan
		echo "mkdir -p ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/reads" >> hic_01_HIChiglassPrepare_single_${CONT_DB}.${slurmID}.plan
		echo "mkdir -p ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/bams" >> hic_01_HIChiglassPrepare_single_${CONT_DB}.${slurmID}.plan
		echo "mkdir -p ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/filter" >> hic_01_HIChiglassPrepare_single_${CONT_DB}.${slurmID}.plan
		echo "mkdir -p ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/matrix" >> hic_01_HIChiglassPrepare_single_${CONT_DB}.${slurmID}.plan		
		echo "mkdir -p ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/ref" >> hic_01_HIChiglassPrepare_single_${CONT_DB}.${slurmID}.plan
		echo "mkdir -p ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/config" >> hic_01_HIChiglassPrepare_single_${CONT_DB}.${slurmID}.plan
		echo "mkdir -p ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/tmp" >> hic_01_HIChiglassPrepare_single_${CONT_DB}.${slurmID}.plan
		
		if [[ -f ${SC_HIC_REF_EXCLUDELIST} ]]
   		then
   			echo "${SEQKIT_PATH} grep -v -f ${SC_HIC_REF_EXCLUDELIST} ${SC_HIC_REF} | sed -e \"s/:/-/g\" ${SC_HIC_REF} > ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/ref/$(basename ${SC_HIC_REF})" >> hic_01_HIChiglassPrepare_single_${CONT_DB}.${slurmID}.plan
   			echo "${SEQKIT_PATH} grep -f ${SC_HIC_REF_EXCLUDELIST} ${SC_HIC_REF} | sed -e \"s/:/-/g\" ${SC_HIC_REF} > ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/ref/exclude.fasta" >> hic_01_HIChiglassPrepare_single_${CONT_DB}.${slurmID}.plan
   		else   			
   			echo "sed -e \"s/:/-/g\" ${SC_HIC_REF} > ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/ref/$(basename ${SC_HIC_REF})" >> hic_01_HIChiglassPrepare_single_${CONT_DB}.${slurmID}.plan
   		fi		
		echo "samtools faidx ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/ref/$(basename ${SC_HIC_REF})" >> hic_01_HIChiglassPrepare_single_${CONT_DB}.${slurmID}.plan
		echo "awk '{print \$1\"\t\"\$2}' ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/ref/$(basename ${SC_HIC_REF}).fai > ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/ref/$(basename ${SC_HIC_REF%.fasta}).chrom.sizes" >> hic_01_HIChiglassPrepare_single_${CONT_DB}.${slurmID}.plan
		echo "bwa index -b 100000000 ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/ref/$(basename ${SC_HIC_REF})" >> hic_01_HIChiglassPrepare_single_${CONT_DB}.${slurmID}.plan
		echo "cp ${configFile} ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/config/$(basename ${configFile%.sh})_$(date '+%Y-%m-%d_%H-%M-%S').sh" >> hic_01_HIChiglassPrepare_single_${CONT_DB}.${slurmID}.plan
		
		echo "samtools $(${CONDA_HIC_ENV} && samtools 2>&1 | grep Version | awk '{print $2}' && conda deactivate)" > hic_01_HIChiglassPrepare_single_${CONT_DB}.${slurmID}.version
		echo "bwa $(${CONDA_HIC_ENV} && bwa 2>&1 | grep Version | awk '{print $2}' && conda deactivate)" >> hic_01_HIChiglassPrepare_single_${CONT_DB}.${slurmID}.version        
	### 02_HiChiglassBwa
	elif [[ ${currentStep} -eq 2 ]]
    then
        ### clean up plans 
        for x in $(ls hic_02_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
                
        if [[ ! -d "${SC_HIC_READS}" ]]
        then
        	(>&2 echo "ERROR - set SC_HIC_READS to directory that contain the PROJECT_ID*.fastq.qz read files")
        	exit 1
   		fi   		   				
        
        if [[ ! -d "${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/reads" ]]
        then
        	(>&2 echo "ERROR - cannot access directory ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/reads!")
        	exit 1
   		fi
   		
   		ref=${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/ref/$(basename ${SC_HIC_REF})
   		
   		if [[ ! -f "${ref}" ]]
        then
        (>&2 echo "ERROR - cannot reference fasta file ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/ref/$(basename ${SC_HIC_REF})!")
        	exit 1
   		fi
   		
   		### link HiC reads into current reads sub directory
   		for x in ${SC_HIC_READS}/${PROJECT_ID}_*_*_R[12].fastq.gz
		do
			if [[ -f ${x} ]]
			then	
				ln -s -r -f ${x} ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/reads 
			fi
		done
   		   				
   		numR1Files=0
		for x in ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/reads/${PROJECT_ID}_*_*_R1.fastq.gz
		do
			if [[ -f ${x} ]]
			then	
				numR1Files=$((${numR1Files}+1))	
			fi
		done
		
		if [[ ${numR1Files} -eq 0 ]]
        then
        	(>&2 echo "ERROR - cannot read HiC R1 files with following pattern: ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/reads/${PROJECT_ID}_*_*_R1.fastq.gz")
        	exit 1
   		fi
   		
   		numR2Files=0
		for x in ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/reads/${PROJECT_ID}_*_*_R2.fastq.gz
		do
			if [[ -f ${x} ]]
			then	
				numR2Files=$((${numR2Files}+1))	
			fi
		done
		
		if [[ ${numR2Files} -eq 0 ]]
        then
        	(>&2 echo "ERROR - cannot read HiC R2 files with following pattern: ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/reads/${PROJECT_ID}_*_*_R2.fastq.gz")
        	exit 1
   		fi
   		
   		if [[ ${numR1Files} -ne ${numR2Files} ]]
        then
        	(>&2 echo "ERROR - HiC R1 files ${numR1Files} does not match R2 files ${numR2Files}")
        	exit 1
   		fi
   		
   		setbwaOptions
   		
		for r1 in ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/reads/${PROJECT_ID}_*_*_R1.fastq.gz
		do
			id=$(dirname ${r1})
			f1=$(basename ${r1})
			f2=$(echo "${f1}" | sed -e "s:_R1.fastq.gz:_R2.fastq.gz:")
			o="${f1%_R1.fastq.gz}"											
			
			if [[ ${numR2Files} -eq 1 ]]
        	then
				echo "bwa mem${CONTIG_BWA_OPT} -SP5M -R \"@RG\tID:${o}\tSM:${PROJECT_ID}_HIC\tLB:${PROJECT_ID}_HIC\tPL:ILLUMINA\tPU:none\" ${ref} ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/reads/${f1} ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/reads/${f2} | samtools view -bhS - > ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/bams/${o}_bwa.bam && samtools view -h ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/bams/${o}_bwa.bam | pairtools parse --nproc-in ${SC_HIC_HIGLASS_PAIRTOOLSTHREADS} --nproc-out ${SC_HIC_HIGLASS_PAIRTOOLSTHREADS} -c ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/ref/$(basename ${SC_HIC_REF%.fasta}).chrom.sizes -o ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/filter/${PROJECT_ID}_allHiC.parsed.pairsam.gz && pairtools sort --tmpdir ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/tmp --memory 8G --nproc ${SC_HIC_HIGLASS_PAIRTOOLSTHREADS} -o ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/filter/${PROJECT_ID}_allHiC.sorted.pairsam.gz ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/filter/${PROJECT_ID}_allHiC.parsed.pairsam.gz"              				
			else
				echo "bwa mem${CONTIG_BWA_OPT} -SP5M -R \"@RG\tID:${o}\tSM:${PROJECT_ID}_HIC\tLB:${PROJECT_ID}_HIC\tPL:ILLUMINA\tPU:none\" ${ref} ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/reads/${f1} ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/reads/${f2} | samtools view -bhS - > ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/bams/${o}_bwa.bam && samtools view -h ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/bams/${o}_bwa.bam | pairtools parse --nproc-in ${SC_HIC_HIGLASS_PAIRTOOLSTHREADS} --nproc-out ${SC_HIC_HIGLASS_PAIRTOOLSTHREADS} -c ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/ref/$(basename ${SC_HIC_REF%.fasta}).chrom.sizes -o ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/filter/${o}_bwa.parsed.pairsam.gz && pairtools sort --tmpdir ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/tmp --memory 8G --nproc ${SC_HIC_HIGLASS_PAIRTOOLSTHREADS} -o ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/filter/${o}_bwa.sorted.pairsam.gz ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/filter/${o}_bwa.parsed.pairsam.gz" 					
			fi			 				 
		done > hic_02_HiChiglassBwa_block_${CONT_DB}.${slurmID}.plan
		
   		echo "bwa $(${CONDA_HIC_ENV} && bwa 2>&1 | grep Version | awk '{print $2}' && conda deactivate)" > hic_02_HiChiglassBwa_block_${CONT_DB}.${slurmID}.version
   		echo "samtools $(${CONDA_HIC_ENV} && samtools 2>&1 | grep Version | awk '{print $2}' && conda deactivate)" >> hic_02_HiChiglassBwa_block_${CONT_DB}.${slurmID}.version                        
	### 03_HiChiglassFilter
	elif [[ ${currentStep} -eq 3 ]]
    then
        ### clean up plans 
        for x in $(ls hic_03_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        if [[ ! -d "${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/bams" ]]
        then
        	(>&2 echo "ERROR - cannot access directory ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/bams!")
        	exit 1
   		fi
   		
   		setPicardOptions
   		   		
   		## if multiple bam files are available (e.g. different Lanes) then merge files prior to markduplicates
   		files=$(ls ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/bams/*_bwa.bam)   		   		   	   		
   		
   		if [[ $(echo $files | wc -w) -eq 1 ]]
   		then        
	        echo "pairtools dedup --nproc-in ${SC_HIC_HIGLASS_PAIRTOOLSTHREADS} --nproc-out ${SC_HIC_HIGLASS_PAIRTOOLSTHREADS} -o ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/filter/${PROJECT_ID}_allHiC.dedup.pairsam.gz ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/filter/${PROJECT_ID}_allHiC.sorted.pairsam.gz"
	    	echo "pairtools select --nproc-in ${SC_HIC_HIGLASS_PAIRTOOLSTHREADS} --nproc-out ${SC_HIC_HIGLASS_PAIRTOOLSTHREADS} '(pair_type == \"UU\") or (pair_type == \"UR\") or (pair_type == \"RU\")' -o ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/filter/${PROJECT_ID}_allHiC.filtered.pairsam.gz ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/filter/${PROJECT_ID}_allHiC.dedup.pairsam.gz"
			echo "pairtools split --nproc-in ${SC_HIC_HIGLASS_PAIRTOOLSTHREADS} --nproc-out ${SC_HIC_HIGLASS_PAIRTOOLSTHREADS} --output-pairs ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/filter/${PROJECT_ID}_allHiC.filtered.pairs.gz --output-sam ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/filter/${PROJECT_ID}_allHiC.filtered.bam ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/filter/${PROJECT_ID}_allHiC.filtered.pairsam.gz"
	    	echo "pairix -f -p pairs ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/filter/${PROJECT_ID}_allHiC.filtered.pairs.gz"    	
					   		   
   		elif [[ $(echo $files | wc -w) -gt 1 ]]
   		then
   			echo "pairtools merge --nproc-in ${SC_HIC_HIGLASS_PAIRTOOLSTHREADS} --nproc-out ${SC_HIC_HIGLASS_PAIRTOOLSTHREADS} -o ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/filter/${PROJECT_ID}_allHiC.sorted.pairsam.gz \$(ls ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/filter/*_bwa.sorted.pairsam.gz)"
	        echo "pairtools dedup --nproc-in ${SC_HIC_HIGLASS_PAIRTOOLSTHREADS} --nproc-out ${SC_HIC_HIGLASS_PAIRTOOLSTHREADS} -o ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/filter/${PROJECT_ID}_allHiC.dedup.pairsam.gz ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/filter/${PROJECT_ID}_allHiC.sorted.pairsam.gz"
	    	echo "pairtools select --nproc-in ${SC_HIC_HIGLASS_PAIRTOOLSTHREADS} --nproc-out ${SC_HIC_HIGLASS_PAIRTOOLSTHREADS} '(pair_type == \"UU\") or (pair_type == \"UR\") or (pair_type == \"RU\")' -o ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/filter/${PROJECT_ID}_allHiC.filtered.pairsam.gz ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/filter/${PROJECT_ID}_allHiC.dedup.pairsam.gz"
			echo "pairtools split --nproc-in ${SC_HIC_HIGLASS_PAIRTOOLSTHREADS} --nproc-out ${SC_HIC_HIGLASS_PAIRTOOLSTHREADS} --output-pairs ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/filter/${PROJECT_ID}_allHiC.filtered.pairs.gz --output-sam ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/filter/${PROJECT_ID}_allHiC.filtered.bam ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/filter/${PROJECT_ID}_allHiC.filtered.pairsam.gz"
	    	echo "pairix -f -p pairs ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/filter/${PROJECT_ID}_allHiC.filtered.pairs.gz"
   		else
   			(>&2 echo "ERROR - cannot find file with following pattern: ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/bams/*_bwaFilt.bam!")
        	exit 1
   	 	fi > hic_03_HiChiglassFilter_single_${CONT_DB}.${slurmID}.plan
        
	### 04_HiChiglassMatrix
	elif [[ ${currentStep} -eq 4 ]]
    then
        ### clean up plans 
        for x in $(ls hic_04_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done        
        
        for i in "${SC_HIC_HIGLASS_COOLERRESOLUTION[@]}"
		do
	      	echo "HDF5_USE_FILE_LOCKING=FALSE cooler cload pairix -p ${SC_HIC_HIGLASS_PAIRTOOLSTHREADS} ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/ref/$(basename ${SC_HIC_REF%.fasta}).chrom.sizes:${i} ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/filter/${PROJECT_ID}_allHiC.filtered.pairs.gz ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/matrix/${PROJECT_ID}_allHiC.output.${i}.cool"   	
    		# normalization - (matrix balancing)
        	echo "HDF5_USE_FILE_LOCKING=FALSE cooler balance ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/matrix/${PROJECT_ID}_allHiC.output.${i}.cool"
        	# aggregation - (for HiGlass view)
        	echo "HDF5_USE_FILE_LOCKING=FALSE cooler zoomify ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/matrix/${PROJECT_ID}_allHiC.output.${i}.cool"
    	done > hic_04_HiChiglassMatrix_single_${CONT_DB}.${slurmID}.plan
	else	
    	(>&2 echo "step ${currentStep} in SC_HIC_TYPE ${SC_HIC_TYPE} not supported")
    	(>&2 echo "valid steps are: ${myTypes[${SC_HIC_TYPE}]}")
    	exit 1
	fi
else
    (>&2 echo "unknown SC_HIC_TYPE ${SC_HIC_TYPE}")
    (>&2 echo "supported types")
    x=0; while [ $x -lt ${#myTypes[*]} ]; do (>&2 echo "${myTypes[${x}]}"); done 
    exit 1
fi

exit 0