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

function setPBalignOptions()
{
    ARROW_PBALIGN_OPT=""
	if [[ -n ${PB_ARROW_PBALIGN_THREADS} && ${PB_ARROW_PBALIGN_THREADS} -gt 0 ]]
	then 
		ARROW_PBALIGN_OPT="${ARROW_PBALIGN_OPT} --nproc ${PB_ARROW_PBALIGN_THREADS}"
	fi	

	if [[ -n ${PB_ARROW_PBALIGN_MINLEN} && ${PB_ARROW_PBALIGN_MINLEN} -gt 0 ]]
	then 
		ARROW_PBALIGN_OPT="${ARROW_PBALIGN_OPT} --minLength ${PB_ARROW_PBALIGN_MINLEN}"
	fi	
	
	### add sa file by default
	if [[ ! -s ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/arrow_in.fasta.sa ]]
	then
		(>&2 echo "missing sa file: ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/arrow_in.fasta.sa")
    	exit 1			
	fi
	ARROW_PBALIGN_OPT="${ARROW_PBALIGN_OPT} --algorithmOptions=\"--sa ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/arrow_in.fasta.sa\""
}

function setArrowOptions()
{
	ARROW_ARROW_OPT=""
	if [[ -n ${PB_ARROW_ARROW_THREADS} && ${PB_ARROW_ARROW_THREADS} -gt 0 ]]
	then 
		ARROW_ARROW_OPT="${ARROW_ARROW_OPT} -j ${PB_ARROW_ARROW_THREADS}"
	fi	
	
	if [[ -n ${PB_ARROW_ARROW_DIPLOID} && ${PB_ARROW_ARROW_DIPLOID} -gt 0 ]]
	then 
		ARROW_ARROW_OPT="${ARROW_ARROW_OPT} --diploid"
	fi
	
	if [[ -n ${PB_ARROW_ARROW_VERBOSE} && ${PB_ARROW_ARROW_VERBOSE} -gt 0 ]]
	then 
		ARROW_ARROW_OPT="${ARROW_ARROW_OPT} -v"
	fi	
	
	if [[ -n ${PB_ARROW_ARROW_REPORTEFFCOV} && ${PB_ARROW_ARROW_REPORTEFFCOV} -gt 0 ]]
	then 
		ARROW_ARROW_OPT="${ARROW_ARROW_OPT} --reportEffectiveCoverage"
	fi
	
	if [[ -n ${PB_ARROW_ARROW_GFFOUT} && ${PB_ARROW_ARROW_GFFOUT} -gt 0 && -n ${PB_ARROW_ARROW_ANNOTATEGFF} && ${PB_ARROW_ARROW_ANNOTATEGFF} -gt 0 ]]
	then 
		ARROW_ARROW_OPT="${ARROW_ARROW_OPT} --annotateGFF"
	fi		
}

myTypes=("1-prepInFasta, 2-pbalign, 3-bamsplit, 4-bamseparate, 5-bamMerge, 6-arrow, 7-statistics")
if [[ ${PB_ARROW_TYPE} -eq 0 ]]
then 
    ### 1-prepInFasta
    if [[ ${currentStep} -eq 1 ]]
    then
        ### clean up plans 
        for x in $(ls arrow_01_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        if [[ ! -f "${PB_ARROW_REFFASTA}" ]]
        then
        	(>&2 echo "ERROR - set PB_ARROW_REFFASTA to input fasta file")
        	exit 1
   		fi
   				
		if [[ ! -d ${PB_ARROW_OUTDIR} ]] 
		then
			(>&2 echo "ERROR - Variable ${PB_ARROW_OUTDIR} is not set or cannot be accessed")
        	exit 1
		fi
		
		echo "if [[ -d ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID} ]]; then mv ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID} ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}_$(date '+%Y-%m-%d_%H-%M-%S'); fi && mkdir ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}" > arrow_01_prepInFasta_single_${CONT_DB}.${slurmID}.plan
		if [[ -n ${PB_ARROW_MAKEUNIQUEHEADER} && ${PB_ARROW_MAKEUNIQUEHEADER} -eq 1 ]]
		then
			echo "cat ${PB_ARROW_REFFASTA} | awk -v count=1 '{if (\$1 ~ /^>/) {print $1\"_\"; count+=1;} else print \$1;}' | sed \"s:|:_:g\" > ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/arrow_in.fasta" >> arrow_01_prepInFasta_single_${CONT_DB}.${slurmID}.plan
		else
			echo "cat ${PB_ARROW_REFFASTA} | awk '{ print \$1}' | sed \"s:|:_:g\" > ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/arrow_in.fasta" >> arrow_01_prepInFasta_single_${CONT_DB}.${slurmID}.plan
		fi
		echo "sawriter ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/arrow_in.fasta" >> arrow_01_prepInFasta_single_${CONT_DB}.${slurmID}.plan
		echo "samtools faidx ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/arrow_in.fasta" >> arrow_01_prepInFasta_single_${CONT_DB}.${slurmID}.plan
		echo "grep -e \">\" ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/arrow_in.fasta | sed -e 's:^>::' > ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/arrow_in.header" >> arrow_01_prepInFasta_single_${CONT_DB}.${slurmID}.plan
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > arrow_01_prepInFasta_single_${CONT_DB}.${slurmID}.version
		echo "samtools $(${CONDA_BASE_ENV} && samtools 2>&1 | grep Version | awk '{print $2}' && conda deactivate)" >> arrow_01_prepInFasta_single_${CONT_DB}.${slurmID}.version
    ### 2-pbalign
    elif [[ ${currentStep} -eq 2 ]]
    then
        ### clean up plans 
        for x in $(ls arrow_02_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        if [[ ! -d ${PB_ARROW_BAM} ]]
        then
        	(>&2 echo "ERROR - Variable ${PB_ARROW_BAM} is not set or cannot be accessed")
        	exit 1
        fi
        
        if [[ ! -f ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/arrow_in.header ]]
        then 
        (>&2 echo "ERROR - missing file ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/arrow_in.header")
	       	exit 1
        fi
        
        # sanity checks
   		numFiles=0 
   		for x in ${PB_ARROW_BAM}/*.subreads.bam   		
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
		
		setPBalignOptions
        
        for x in ${PB_ARROW_BAM}/*.subreads.bam   		
   		do
        	name=$(basename ${x%.subreads.bam})
        	
        	logfile=""
    		if [[ -n ${PB_ARROW_PBALIGN_LOGFILE} && ${PB_ARROW_PBALIGN_LOGFILE} -ne 0 ]]
    		then
    			logfile=" --log-file ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${name}.pbalign.log"	
    		fi
    		unalignFile=""
    		if [[ -n ${PB_ARROW_PBALIGN_UNALIGNFILE} && ${PB_ARROW_PBALIGN_UNALIGNFILE} -ne 0 ]]
    		then
    			unalignFile=" --unaligned ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${name}.pbalignUnalign.txt"	
    		fi
    		
    		echo "pbalign${logfile}${unalignFile}${ARROW_PBALIGN_OPT} ${x} ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/arrow_in.fasta ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${name}.pbalign.bam"        	
    	done > arrow_02_pbalign_block_${CONT_DB}.${slurmID}.plan 
    	echo "pbalign $(${CONDA_BASE_ENV} && pbalign --version && conda deactivate)" > arrow_02_pbalign_block_${CONT_DB}.${slurmID}.version
		echo "samtools $(${CONDA_BASE_ENV} && samtools 2>&1 | grep Version | awk '{print $2}' && conda deactivate)" >> arrow_02_pbalign_block_${CONT_DB}.${slurmID}.version
		echo "blasr $(${CONDA_BASE_ENV} && blasr --version | awk '{print $2}' && conda deactivate)">> arrow_02_pbalign_block_${CONT_DB}.${slurmID}.version
    ### 3-bamsplit
    elif [[ ${currentStep} -eq 3 ]]
    then
        ### clean up plans 
        for x in $(ls arrow_03_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        if [[ ! -d ${PB_ARROW_BAM} ]]
        then
        	(>&2 echo "ERROR - Variable ${PB_ARROW_BAM} is not set or cannot be accessed")
        	exit 1
        fi
        
        if [[ ! -f ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/arrow_in.header ]]
        then 
        	(>&2 echo "ERROR - missing file ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/arrow_in.header")
	       	exit 1
        fi
        
        # sanity checks
   		numFiles=0 
		for x in ${PB_ARROW_BAM}/*.subreads.bam   		
   		do
        	name=$(basename ${x%.subreads.bam})
        	file=${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${name}.pbalign.bam
        	if [[ ! -f ${file} || ! -s ${file} ]]
   			then
   				(>&2 echo "WARNING - file ${file} not available or empty")
   			else
   				numFiles=$((${numFiles}+1))
   			fi      						
   		done

		if [[ ${numFiles} -eq 0 ]]
		then
			(>&2 echo "ERROR - no input bam file found")
	       	exit 1	
		fi
		
		for x in ${PB_ARROW_BAM}/*.subreads.bam   		
   		do
        	name=$(basename ${x%.subreads.bam})
        	file=${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${name}.pbalign.bam
        	echo "bamtools split -in  ${file} -reference"			   			   		
		done > arrow_03_bamtools_block_${CONT_DB}.${slurmID}.plan
		echo "$(${CONDA_BASE_ENV} && bamtools --version | grep bamtools && conda deactivate)" > arrow_03_bamtools_block_${CONT_DB}.${slurmID}.version
    ### 4-bamseparate
    elif [[ ${currentStep} -eq 4 ]]
    then
        ### clean up plans 
        for x in $(ls arrow_04_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
	    if [[ ! -f ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/arrow_in.header ]]
	    then 
	    	(>&2 echo "ERROR - missing file ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/arrow_in.header")
	       	exit 1
	    fi
	    
		echo "for x in \$(cat ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/arrow_in.header); do mkdir -p ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/\$x; mv ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/*.REF_\${x}.bam ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/\${x}; find \$(pwd)/${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/\${x} -name \"*.bam\" > ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/\${x}/in.fof; done" > arrow_04_bamseparate_single_${CONT_DB}.${slurmID}.plan	    
	### 5-bamMerge 
    elif [[ ${currentStep} -eq 5 ]]
    then
        ### clean up plans 
        for x in $(ls arrow_05_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
                
	    if [[ ! -f ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/arrow_in.header ]]
	    then 
	    	(>&2 echo "ERROR - missing file ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/arrow_in.header")
	       	exit 1
	    fi
	    
		# sanity checks
   		numFiles=0 
   		for x in $(cat ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/arrow_in.header)   		
   		do
   			if [[ ! -f ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${x}/in.fof || ! -s ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${x}/in.fof ]]
   			then
   				(>&2 echo "WARNING - file ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${x}/in.fof not available or empty")
   			else
   				numFiles=$((${numFiles}+1))
   			fi      						
   		done

		if [[ ${numFiles} -eq 0 ]]
		then
			(>&2 echo "ERROR - no splitted reference based alignment bam files found")
	       	exit 1	
		fi
		
		for x in $(cat ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/arrow_in.header)   		
   		do
   			echo "bamtools merge -list ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${x}/in.fof -out ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${x}/ALL_${x}.bam && if [[ -s ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${x}/in.fof ]]; then xargs rm < ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${x}/in.fof; fi"
		done > arrow_05_bamtools_block_${CONT_DB}.${slurmID}.plan
		echo "$(${CONDA_BASE_ENV} && bamtools --version | grep bamtools && conda deactivate)" >arrow_05_bamtools_block_${CONT_DB}.${slurmID}.version
	### 6-arrow 
    elif [[ ${currentStep} -eq 6 ]]
    then
        ### clean up plans 
        for x in $(ls arrow_06_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        if [[ ! -f ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/arrow_in.header ]]
	    then 
	    	(>&2 echo "ERROR - missing file ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/arrow_in.header")
	       	exit 1
	    fi
	    
		# sanity checks
   		numFiles=0 
   		for x in $(cat ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/arrow_in.header)   		
   		do
   			if [[ ! -f ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${x}/ALL_${x}.bam || ! -s ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${x}/ALL_${x}.bam ]]
   			then
   				(>&2 echo "WARNING - file ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${x}/ALL_${x}.bam not available or empty")
   			else
   				numFiles=$((${numFiles}+1))
   			fi      						
   		done

		if [[ ${numFiles} -eq 0 ]]
		then
			(>&2 echo "ERROR - no merged reference based alignment bam file found")
	       	exit 1	
		fi
		
		setArrowOptions
		
		for x in $(cat ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/arrow_in.header)   		
   		do
   			gff=""
   			if [[ -n ${PB_ARROW_ARROW_GFFOUT} && ${PB_ARROW_ARROW_GFFOUT} -ne 0 ]]
   			then
   				gff=" -o ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${x}/ALL_${x}.arrow.gff"	
   			fi 
   			vcf=""
   			if [[ -n ${PB_ARROW_ARROW_VCFOUT} && ${PB_ARROW_ARROW_VCFOUT} -ne 0 ]]
   			then
   				vcf=" -o ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${x}/ALL_${x}.arrow.vcf"	
   			fi
   			fq=""
   			if [[ -n ${PB_ARROW_ARROW_FQOUT} && ${PB_ARROW_ARROW_FQOUT} -ne 0 ]]
   			then
   				fq=" -o ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${x}/ALL_${x}.arrow.fq"	
   			fi
   			
   			# check if bam is empty !!!!!
   			if [[ $(samtools view -H ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${x}/ALL_${x}.bam | wc -l) -eq 0 ]]
   			then
   				(>&2 echo "bam file ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${x}/ALL_${x}.bam is EMPTY. Skip it!")
   			else
   				echo "bamtools index -in ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${x}/ALL_${x}.bam && pbindex ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${x}/ALL_${x}.bam && arrow${ARROW_ARROW_OPT}${gff}${vcf}${fq} -r ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/arrow_in.fasta -w ${x} -o ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${x}/ALL_${x}.arrow.fa ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${x}/ALL_${x}.bam"	
   			fi
		done > arrow_06_arrow_block_${CONT_DB}.${slurmID}.plan	
		echo "$(${CONDA_BASE_ENV} && bamtools --version | grep bamtools && conda deactivate)" > arrow_06_arrow_block_${CONT_DB}.${slurmID}.version
		echo "pbindex $(${CONDA_BASE_ENV} && pbindex --version && conda deactivate)" >> arrow_06_arrow_block_${CONT_DB}.${slurmID}.version
		echo "arrow $(${CONDA_BASE_ENV} && arrow --version && conda deactivate)" >> arrow_06_arrow_block_${CONT_DB}.${slurmID}.version				
	### 7-statistics 
    elif [[ ${currentStep} -eq 7 ]]
    then
        ### clean up plans 
        for x in $(ls arrow_07_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
                	
        if [[ -n ${PB_ARROW_FULLSTATS} && ${PB_ARROW_FULLSTATS} -gt 0 ]]
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
    	echo "${SUBMIT_SCRIPTS_PATH}/assemblyStats.sh ${configFile} 9" > arrow_07_statistics_single_${CONT_DB}.${slurmID}.plan
    	git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD > arrow_07_statistics_single_${CONT_DB}.${slurmID}.version
    else
        (>&2 echo "step ${currentStep} in PB_ARROW_TYPE ${PB_ARROW_TYPE} not supported")
        (>&2 echo "valid steps are: ${myTypes[${PB_ARROW_TYPE}]}")
        exit 1            
    fi    
else
    (>&2 echo "unknown PB_ARROW_TYPE ${PB_ARROW_TYPE}")
    (>&2 echo "supported types")
    x=0; while [ $x -lt ${#myTypes[*]} ]; do (>&2 echo "${myTypes[${x}]}"); done 
    exit 1
fi

exit 0