#!/bin/bash -e

configFile=$1
currentStep=$2
slurmID=$3

if [[ ! -f ${configFile} ]]
then 
    (>&2 echo "cannot access config file ${configFile}")
    exit 1
fi

function getNumOfDbBlocks()
{
    db=$1
    if [[ ! -f $db ]]
    then
        (>&2 echo "database $db not found")
        exit 1
    fi

    blocks=$(grep block $db | awk '{print $3}')
    if [[ ! -n $blocks ]]
    then 
        (>&2 echo "database $db has not been partitioned. Run DBsplit first!")
        exit 1
    fi 
    echo ${blocks}
}

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
myCWD=$(pwd)
myTypes=("1-prepInFasta, 2-createMinimap2RefIndex, 3-minimap2, 4-bamMerge, 5-readCovHist, 6-contigCovHist, 7-purgeHaplotigs, 8-statistics",
"01_PDprepInput, 02_PDminimap2, 03_PDcalcuts, 04_PDminimap2, 05_purgedups, 06_statistics", "01_TCPrepInput, 02_TCDBdust, 03_TCdatander, 04_TCCatrack, 05_TCdaligner, 06_TCLAmerge, 07_TCLAfilterChain, 08_TCLAmerge, 09_TCCTtrim, 10_TCstatistics")
if [[ ${CT_PURGEHAPLOTIGS_TYPE} -eq 0 ]]
then 
    ### 1-prepInFasta
    if [[ ${currentStep} -eq 1 ]]
    then
        ### clean up plans 
        for x in $(ls purgeHaplotigs_01_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        if [[ ! -f "${CT_PURGEHAPLOTIGS_INFASTA}" ]]
        then
        	(>&2 echo "ERROR - set CT_PURGEHAPLOTIGS_INFASTA to input fasta file")
        	exit 1
   		fi
   		
   		echo "if [[ -d ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID} ]]; then mv ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID} ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}_$(date '+%Y-%m-%d_%H-%M-%S'); fi && mkdir ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}" > purgeHaplotigs_01_prepInFasta_single_${CONT_DB}.${slurmID}.plan
		echo "ln -s -r ${CT_PURGEHAPLOTIGS_INFASTA} ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}" >> purgeHaplotigs_01_prepInFasta_single_${CONT_DB}.${slurmID}.plan
	### 2-createMinimap2RefIndex
    elif [[ ${currentStep} -eq 2 ]]
    then
    	### clean up plans 
        for x in $(ls purgeHaplotigs_02_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        ref=$(basename ${CT_PURGEHAPLOTIGS_INFASTA%.fasta})
        
        echo "minimap2 -t ${CT_PURGEHAPLOTIGS_MINIMAP2IDXTHREADS} -x map-pb -d ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}/${ref}.idx ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}/${ref}.fasta" > purgeHaplotigs_02_createMinimap2RefIndex_single_${CONT_DB}.${slurmID}.plan
        echo "minimap2 $(${PURGEHAPLOTIGS_ENV} && minimap2 --version && ${PURGEHAPLOTIGS_ENV_DEACT})" > purgeHaplotigs_02_createMinimap2RefIndex_single_${CONT_DB}.${slurmID}.version
		echo "samtools $(${PURGEHAPLOTIGS_ENV} && samtools 2>&1 | grep Version | awk '{print $2}' && ${PURGEHAPLOTIGS_ENV_DEACT})" >> purgeHaplotigs_02_createMinimap2RefIndex_single_${CONT_DB}.${slurmID}.version
    ### 3-minimap2
    elif [[ ${currentStep} -eq 3 ]]
    then
        ### clean up plans 
        for x in $(ls purgeHaplotigs_03_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
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
   		for x in ${CT_PURGEHAPLOTIGS_PACBIOFASTA}/*.subreads.fa.gz   		
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
		
		if [[ ! -f ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}/${ref}.idx ]]
		then 
			(>&2 echo "ERROR - file ${ref}.idx not available. Create an reference index first!")
			exit 1
		fi
		
        for x in ${CT_PURGEHAPLOTIGS_PACBIOFASTA}/*.subreads.fa.gz   		
   		do
        	name=$(basename ${x%.subreads.fa.gz})
        	    		
    	echo "minimap2 -a -x map-pb -t ${CT_PURGEHAPLOTIGS_MINIMAP2ALNTHREADS} ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}/${ref}.idx ${x} | samtools view -h - | samtools sort -@ ${CT_PURGEHAPLOTIGS_SAMTOOLSTHREADS} -m ${CT_PURGEHAPLOTIGS_SAMTOOLSMEM}G -o ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}/${ref}_${name}_minimap2.sort.bam -T /tmp/${ref}_${name}_minimap2.sort.tmp && samtools index -@ ${CT_PURGEHAPLOTIGS_SAMTOOLSTHREADS} ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}/${ref}_${name}_minimap2.sort.bam"        	
		done > purgeHaplotigs_03_minimap2_block_${CONT_DB}.${slurmID}.plan 
    	echo "minimap2 $(${PURGEHAPLOTIGS_ENV} && minimap2 --version && ${PURGEHAPLOTIGS_ENV_DEACT})" > purgeHaplotigs_03_minimap2_block_${CONT_DB}.${slurmID}.version
		echo "samtools $(${PURGEHAPLOTIGS_ENV} && samtools 2>&1 | grep Version | awk '{print $2}' && ${PURGEHAPLOTIGS_ENV_DEACT})" >> purgeHaplotigs_03_minimap2_block_${CONT_DB}.${slurmID}.version
	### 4-bamMerge
    elif [[ ${currentStep} -eq 4 ]]
    then
        ### clean up plans 
        for x in $(ls purgeHaplotigs_04_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        if [[ ! -d ${CT_PURGEHAPLOTIGS_PACBIOFASTA} ]]
        then
        	(>&2 echo "ERROR - Variable ${CT_PURGEHAPLOTIGS_PACBIOFASTA} is not set or cannot be accessed")
        	exit 1
        fi
        
        ref=$(basename ${CT_PURGEHAPLOTIGS_INFASTA%.fasta})
		
		if [[ ! -f ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}/${ref}.idx ]]
		then 
			(>&2 echo "ERROR - file ${ref}.idx not available. Create an reference index first!")
			exit 1
		fi
        
                
    	# create input bam file list 
   		for x in ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}/${ref}_*_minimap2.sort.bam   		
   		do
   			if [[ ! -f ${x} || ! -s ${x} ]]
   			then
   				(>&2 echo "WARNING - file ${x} not available or empty")
   			else
   				echo "${x}" 
   			fi      						
		done > ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}/bamlist.fofn
        	    		
    	echo "samtools merge -b ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}/bamlist.fofn -@ ${CT_PURGEHAPLOTIGS_SAMTOOLSTHREADS} ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}/${ref}_minimap2.sort.bam" > purgeHaplotigs_04_bamMerge_single_${CONT_DB}.${slurmID}.plan
    	echo "samtools index -@ ${CT_PURGEHAPLOTIGS_SAMTOOLSTHREADS} ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}/${ref}_minimap2.sort.bam" >> purgeHaplotigs_04_bamMerge_single_${CONT_DB}.${slurmID}.plan
   		echo "samtools $(${PURGEHAPLOTIGS_ENV} && samtools 2>&1 | grep Version | awk '{print $2}' && ${PURGEHAPLOTIGS_ENV_DEACT})" > purgeHaplotigs_04_bamMerge_single_${CONT_DB}.${slurmID}.version
   	### 5-readCovHist
    elif [[ ${currentStep} -eq 5 ]]
    then
        ### clean up plans 
        for x in $(ls purgeHaplotigs_05_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        if [[ ! -d ${CT_PURGEHAPLOTIGS_PACBIOFASTA} ]]
        then
        	(>&2 echo "ERROR - Variable ${CT_PURGEHAPLOTIGS_PACBIOFASTA} is not set or cannot be accessed")
        	exit 1
        fi
        
        ref=$(basename ${CT_PURGEHAPLOTIGS_INFASTA%.fasta})
		
		if [[ ! -f ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}/${ref}.idx ]]
		then 
			(>&2 echo "ERROR - file ${ref}.idx not available. Create an reference index first!")
			exit 1
		fi
        
		if [[ ! -f ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}/${ref}_minimap2.sort.bam ]] 
        then 
        	(>&2 echo "ERROR - file ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}/${ref}_minimap2.sort.bam not available.")
			exit 1
    	fi
    	
    	echo "purge_haplotigs readhist -b ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}/${ref}_minimap2.sort.bam -g ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}/${ref}.fasta -t ${CT_PURGEHAPLOTIGS_THREADS}" > purgeHaplotigs_05_readCovHist_single_${CONT_DB}.${slurmID}.plan
    	echo "purge_haplotigs $(${PURGEHAPLOTIGS_ENV} && conda list purge_haplotigs | grep -e "^purge_haplotigs" | awk '{print $2}' && ${PURGEHAPLOTIGS_ENV_DEACT})" > purgeHaplotigs_05_readCovHist_single_${CONT_DB}.${slurmID}.version	   						
    else
        (>&2 echo "step ${currentStep} in CT_PURGEHAPLOTIGS_TYPE ${CT_PURGEHAPLOTIGS_TYPE} not supported")
        (>&2 echo "valid steps are: ${myTypes[${CT_PURGEHAPLOTIGS_TYPE}]}")
        exit 1            
    fi    		
elif [[ ${CT_PURGEHAPLOTIGS_TYPE} -eq 1 ]]
then 
	if [[ -z ${PBBIOCONDA_ENV} ]]
	then 
		(>&2 echo "ERROR - var PBBIOCONDA_ENV needs to be set! We need a conda env with folling tools: bam2fasta")
		exit 1
	fi  

	if [[ -z ${PURGEDUPS_PATH} || ! -f ${PURGEDUPS_PATH}/bin/purge_dups ]]
	then 
		(>&2 echo "ERROR - var ${PURGEDUPS_PATH} needs to be set! We need a installation dir like: \${PURGEDUPS_PATH}/bin/purge_dups!")
		exit 1
	fi  

    ### 01_PDprepInput			("01_PDprepInput, 02_PDminimap2, 03_PDcalcuts, 04_PDminimap2, 05_purgedups, 06_statistics")
    if [[ ${currentStep} -eq 1 ]]		
    then
        ### clean up plans 
        for x in $(ls purgeHaplotigs_01_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        if [[ ! -f "${CT_PURGEHAPLOTIGS_INFASTA}" ]]
        then
        	(>&2 echo "ERROR - set CT_PURGEHAPLOTIGS_INFASTA to input fasta file")
        	exit 1
   		fi
   		
   		echo "if [[ -d ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID} ]]; then mv ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID} ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}_$(date '+%Y-%m-%d_%H-%M-%S'); fi && mkdir ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}" > purgeHaplotigs_01_prepInFasta_single_${CONT_DB}.${slurmID}.plan
		echo "ln -s -r ${CT_PURGEHAPLOTIGS_INFASTA} ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}" >> purgeHaplotigs_01_prepInFasta_single_${CONT_DB}.${slurmID}.plan
    ### 02_PDminimap2			("01_PDprepInput, 02_PDminimap2, 03_PDcalcuts, 04_PDminimap2, 05_purgedups, 06_statistics")
    elif [[ ${currentStep} -eq 2 ]]
    then
        ### clean up plans 
        for x in $(ls purgeHaplotigs_02_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done

		if [[ ! -f ${CT_PURGEHAPLOTIGS_INFASTA} ]]
		then 
			(>&2 echo "ERROR - No reference fasta present. Var CT_PURGEHAPLOTIGS_INFASTA must be set!")
		fi
        
        if [[ -d ${CT_PURGEHAPLOTIGS_PACBIOFASTA} ]]
        then
			# sanity checks
			numFiles=0 
			for x in ${CT_PURGEHAPLOTIGS_PACBIOFASTA}/*.subreads.fasta
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
				(>&2 echo "ERROR - no input *.subreads.fasta files found")
				exit 1	
			fi
			
			ref=$(basename ${CT_PURGEHAPLOTIGS_INFASTA%.fasta})
					
			for x in ${CT_PURGEHAPLOTIGS_PACBIOFASTA}/*.subreads.fasta	
			do
				name=$(basename ${x%.subreads.fasta})
				echo "minimap2 -x map-pb -t ${CT_PURGEHAPLOTIGS_MINIMAP2ALNTHREADS} ${CT_PURGEHAPLOTIGS_INFASTA} ${x} | gzip -c - > ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}/${ref}_${name}_minimap2.sort.paf.gz"        	
			done > purgeHaplotigs_02_PDminimap2_block_${CONT_DB}.${slurmID}.plan 
		else 
        	(>&2 echo "ERROR - Variable ${CT_PURGEHAPLOTIGS_PACBIOFASTA} is not set. Extract subreads on the fly.")
			
			numFiles=0 
			for x in ${DB_PATH}/*.subreads.bam
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
				(>&2 echo "ERROR - no input *.subreads.bam files found")
				exit 1	
			fi

			ref=$(basename ${CT_PURGEHAPLOTIGS_INFASTA%.fasta})
					
			for x in ${DB_PATH}/*.subreads.bam	
			do
				name=$(basename ${x%.subreads.bam})
				echo -n "${PBBIOCONDA_ENV} && bam2fasta -o ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}/${name%.bam}.fa.gz ${x} "
				echo -e "minimap2 -x map-pb -t ${CT_PURGEHAPLOTIGS_MINIMAP2ALNTHREADS} ${CT_PURGEHAPLOTIGS_INFASTA} ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}/${name%.bam}.fa.gz | gzip -c - > ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}/${ref}_${name}_minimap2.sort.paf.gz && conda deactivate"        	
			done > purgeHaplotigs_02_PDminimap2_block_${CONT_DB}.${slurmID}.plan 
        fi
       	echo "minimap2 $(${PBBIOCONDA_ENV}) && minimap2 --version && conda deactivate" > purgeHaplotigs_02_PDminimap2_block_${CONT_DB}.${slurmID}.version	
		echo "$(${PBBIOCONDA_ENV}) && bam2fasta --version && conda deactivate" > purgeHaplotigs_02_PDminimap2_block_${CONT_DB}.${slurmID}.version
	### 03_PDcalcuts				("01_PDprepInput, 02_PDminimap2, 03_PDcalcuts, 04_PDminimap2, 05_purgedups, 06_statistics")
    elif [[ ${currentStep} -eq 3 ]]
    then
        ### clean up plans 
        for x in $(ls purgeHaplotigs_03_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done

		ref=$(basename ${CT_PURGEHAPLOTIGS_INFASTA%.fasta})

		echo "${PURGEDUPS_PATH}/bin/pbcstat -O ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID} ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}/${ref}_*_minimap2.sort.paf.gz" > purgeHaplotigs_03_PDcalcuts_single_${CONT_DB}.${slurmID}.plan 
		echo "${PURGEDUPS_PATH}/bin/calcuts ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}/PB.stat > ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}/cutoffs 2> ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}/calcults.log" >> purgeHaplotigs_03_PDcalcuts_single_${CONT_DB}.${slurmID}.plan 

		echo "purge_dups $(${PURGEDUPS_PATH}/bin/purge_dups -h 2>&1 | grep Version)" > purgeHaplotigs_03_PDcalcuts_single_${CONT_DB}.${slurmID}.version   	
   	### 04_PDminimap2 				("01_PDprepInput, 02_PDminimap2, 03_PDcalcuts, 04_PDminimap2, 05_purgedups, 06_statistics")
    elif [[ ${currentStep} -eq 4  ]]
    then
        ### clean up plans 
        for x in $(ls purgeHaplotigs_04_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
                
        ref=$(basename ${CT_PURGEHAPLOTIGS_INFASTA%.fasta})
		
		addOpt=""
		if [[ -n ${CT_PURGEHAPLOTIGS_MINIMAP2_PRESET} ]]
		then 
			addOpt="${addOpt} -x ${CT_PURGEHAPLOTIGS_MINIMAP2_PRESET}"
		else
			addOpt="${addOpt} -x asm10"
		fi

		echo "${PURGEDUPS_PATH}/bin/split_fa ${CT_PURGEHAPLOTIGS_INFASTA} > ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}/${ref}_split.fasta" > purgeHaplotigs_04_PDminimap2_single_${CONT_DB}.${slurmID}.plan 
		echo "${PBBIOCONDA_ENV} && minimap2${addOpt} -t ${CT_PURGEHAPLOTIGS_MINIMAP2ALNTHREADS} -DP ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}/${ref}_split.fasta ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}/${ref}_split.fasta | gzip -c - > ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}/${ref}_split.self.paf.gz && conda deactivate" >> purgeHaplotigs_04_PDminimap2_single_${CONT_DB}.${slurmID}.plan 

		echo "purge_dups (split_fa) $(${PURGEDUPS_PATH}/bin/purge_dups -h 2>&1 | grep Version)" > purgeHaplotigs_04_PDminimap2_single_${CONT_DB}.${slurmID}.version 
    	echo "minimap2 $(${PBBIOCONDA_ENV}) && minimap2 --version && conda deactivate" >> purgeHaplotigs_04_PDminimap2_single_${CONT_DB}.${slurmID}.version 
   	### 05_purgedups 				("01_PDprepInput, 02_PDminimap2, 03_PDcalcuts, 04_PDminimap2, 05_purgedups, 06_statistics")
    elif [[ ${currentStep} -eq 5  ]]
    then
        ### clean up plans 
        for x in $(ls purgeHaplotigs_05_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
                
        ref=$(basename ${CT_PURGEHAPLOTIGS_INFASTA%.fasta})
		
		echo "${PURGEDUPS_PATH}/bin/purge_dups -2 -T ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}/cutoffs -c ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}/PB.base.cov ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}/${ref}_split.self.paf.gz > ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}/dups.bed 2> ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}/purge_dups.log" > purgeHaplotigs_05_purgedups_single_${CONT_DB}.${slurmID}.plan 
		echo "${PURGEDUPS_PATH}/bin/get_seqs -p ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}/${ref%.fasta} ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}/dups.bed ${CT_PURGEHAPLOTIGS_INFASTA}" >> purgeHaplotigs_05_purgedups_single_${CONT_DB}.${slurmID}.plan 

		echo "purge_dups (get_seqs) $(${PURGEDUPS_PATH}/bin/purge_dups -h 2>&1 | grep Version)" > purgeHaplotigs_05_purgedups_single_${CONT_DB}.${slurmID}.version  
		echo "purge_dups (purge_dups) $(${PURGEDUPS_PATH}/bin/purge_dups -h 2>&1 | grep Version)" >> purgeHaplotigs_05_purgedups_single_${CONT_DB}.${slurmID}.version  
	elif [[ ${currentStep} -eq 6  ]]
    then
        ### clean up plans 
        for x in $(ls purgeHaplotigs_06_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
		
		echo still todo 
    else
        (>&2 echo "step ${currentStep} in CT_PURGEHAPLOTIGS_TYPE ${CT_PURGEHAPLOTIGS_TYPE} not supported")
        (>&2 echo "valid steps are: ${myTypes[${CT_PURGEHAPLOTIGS_TYPE}]}")
        exit 1            
    fi    	
elif [[ ${CT_PURGEHAPLOTIGS_TYPE} -eq 2 ]]
then 
    ### 01_TCPrepInput			("01_TCPrepInput, 02_TCDBdust, 03_TCdatander, 04_TCCatrack, 05_TCdaligner, 06_TCLAmerge, 07_TCLAfilterChain, 08_TCLAmerge, 09_TCCTtrim, 10_TCstatistics")
    if [[ ${currentStep} -eq 1 ]]		
    then
		### clean up plans 
        for x in $(ls purgeHaplotigs_01_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done

		ref=$(basename ${CT_PURGEHAPLOTIGS_INFASTA%.fasta})

		echo "if [[ -d ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID} ]]; then mv ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID} ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}_$(date '+%Y-%m-%d_%H-%M-%S'); fi && mkdir ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}" > purgeHaplotigs_01_TCPrepInput_single_${CONT_DB}.${slurmID}.plan
		echo "ln -s -r ${CT_PURGEHAPLOTIGS_INFASTA} ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}" >> purgeHaplotigs_01_TCPrepInput_single_${CONT_DB}.${slurmID}.plan
		echo "seqkit split -i ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}/${ref}.fasta" >> purgeHaplotigs_01_TCPrepInput_single_${CONT_DB}.${slurmID}.plan
		echo "grep -e \">\" ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}/${ref}.fasta | tr -d \">\" > ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}/${ref}.header" >> purgeHaplotigs_01_TCPrepInput_single_${CONT_DB}.${slurmID}.plan
		echo "for x in \$(cat ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}/${ref}.header); do len=\$(seqkit stats -T --quiet ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}/${ref}.fasta.split/${ref}.id_\${x}.fasta | awk '{if(NR==2)print \$5}'); awk -v l=\${len} '{if (substr(\$1,1,1) ~ /^>/ ) {print \$1\"/1/0_\"l} else print \$0}' ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}/${ref}.fasta.split/${ref}.id_\${x}.fasta > ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}/${ref}.fasta.split/\${x}.fasta; done" >> purgeHaplotigs_01_TCPrepInput_single_${CONT_DB}.${slurmID}.plan

		echo "cd ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID} && for x in \$(cat ${ref}.header); do ${MARVEL_PATH}/bin/FA2db -v -x0 ${PROJECT_ID}_CT_M ${ref}.fasta.split/\${x}.fasta; done && cd ${myCWD}" >> purgeHaplotigs_01_TCPrepInput_single_${CONT_DB}.${slurmID}.plan
		echo "cd ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID} && for x in \$(cat ${ref}.header); do ${DAZZLER_PATH}/bin/fasta2DB -v ${PROJECT_ID}_CT_Z ${ref}.fasta.split/\${x}.fasta; done && cd ${myCWD}" >> purgeHaplotigs_01_TCPrepInput_single_${CONT_DB}.${slurmID}.plan
		echo "cd ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID} && ${MARVEL_PATH}/bin/DBsplit -s50 ${PROJECT_ID}_CT_M && cd ${myCWD}"  >> purgeHaplotigs_01_TCPrepInput_single_${CONT_DB}.${slurmID}.plan
		echo "cd ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID} && ${DAZZLER_PATH}/bin/DBsplit -s50 ${PROJECT_ID}_CT_Z && cd ${myCWD}" >> purgeHaplotigs_01_TCPrepInput_single_${CONT_DB}.${slurmID}.plan

		echo "$(seqkit version)" > purgeHaplotigs_01_TCPrepInput_single_${CONT_DB}.${slurmID}.version
		echo $(awk --version | head -n 1) >> purgeHaplotigs_01_TCPrepInput_single_${CONT_DB}.${slurmID}.version
    	echo "DAmar $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" >> purgeHaplotigs_01_TCPrepInput_single_${CONT_DB}.${slurmID}.version
    	echo "DAZZLER $(git --git-dir=${DAZZLER_SOURCE_PATH}/DAZZ_DB/.git rev-parse --short HEAD)" >> purgeHaplotigs_01_TCPrepInput_single_${CONT_DB}.${slurmID}.version
	### 02_TCDBdust			("01_TCPrepInput, 02_TCDBdust, 03_TCdatander, 04_TCCatrack, 05_TCdaligner, 06_TCLAmerge, 07_TCLAfilterChain, 08_TCLAmerge, 09_TCCTtrim, 10_TCstatistics")
    elif [[ ${currentStep} -eq 2 ]]		
    then
		### clean up plans 
        for x in $(ls purgeHaplotigs_02_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
		
		### find and set DBdust options 
		##TODO setDBdustOptions
        ### create DBdust commands 
		nblocks=$(getNumOfDbBlocks ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}/${PROJECT_ID}_CT_M.db)
        for x in $(seq 1 ${nblocks})
        do 
           	echo "cd ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID} && ${MARVEL_PATH}/bin/DBdust ${PROJECT_ID}_CT_M.${x} && ${DAZZLER_PATH}/bin/DBdust ${PROJECT_ID}_CT_Z.${x} && cd ${myCWD}"
    	done > purgeHaplotigs_02_TCDBdust_block_${CONT_DB}.${slurmID}.plan
        echo "DAmar $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > purgeHaplotigs_02_TCDBdust_block_${CONT_DB}.${slurmID}.version
    	echo "DAZZLER $(git --git-dir=${DAZZLER_SOURCE_PATH}/DAZZ_DB/.git rev-parse --short HEAD)" >> purgeHaplotigs_02_TCDBdust_block_s${CONT_DB}.${slurmID}.version        	
	### 03_TCdatander			("01_TCPrepInput, 02_TCDBdust, 03_TCdatander, 04_TCCatrack, 05_TCdaligner, 06_TCLAmerge, 07_TCLAfilterChain, 08_TCLAmerge, 09_TCCTtrim, 10_TCstatistics")
    elif [[ ${currentStep} -eq 3 ]]		
    then
		### clean up plans 
        for x in $(ls purgeHaplotigs_03_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done

		 ### find and set datander options 
    	CT_PURGEHAPLOTIGS_DATANDER_OPT=""
    	if [[ -n ${CT_PURGEHAPLOTIGS_DATANDER_THREADS} ]]
    	then
        	CT_PURGEHAPLOTIGS_DATANDER_OPT="${CT_PURGEHAPLOTIGS_DATANDER_OPT} -T${CT_PURGEHAPLOTIGS_DATANDER_THREADS}"
    	fi
    	if [[ -n ${CT_PURGEHAPLOTIGS_DATANDER_MINLEN} ]]
    	then
        	CT_PURGEHAPLOTIGS_DATANDER_OPT="${CT_PURGEHAPLOTIGS_DATANDER_OPT} -l${CT_PURGEHAPLOTIGS_DATANDER_MINLEN}"
    	fi

		### create datander commands
		nblocks=$(getNumOfDbBlocks ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}/${PROJECT_ID}_CT_M.db)
        for x in $(seq 1 ${nblocks})
        do 
            echo "cd ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID} && PATH=${DAZZLER_PATH}/bin:\${PATH} ${DAZZLER_PATH}/bin/datander${CT_PURGEHAPLOTIGS_DATANDER_OPT} ${PROJECT_ID}_CT_Z.${x} && cd ${myCWD}"
		done > purgeHaplotigs_03_TCdatander_block_${CONT_DB}.${slurmID}.plan
        echo "DAZZLER datander $(git --git-dir=${DAZZLER_SOURCE_PATH}/DAMASKER/.git rev-parse --short HEAD)" > purgeHaplotigs_03_TCdatander_block_${CONT_DB}.${slurmID}.version
 	### 04_TCCatrack			("01_TCPrepInput, 02_TCDBdust, 03_TCdatander, 04_TCCatrack, 05_TCdaligner, 06_TCLAmerge, 07_TCLAfilterChain, 08_TCLAmerge, 09_TCCTtrim, 10_TCstatistics")
    elif [[ ${currentStep} -eq 4 ]]		
    then
		### clean up plans 
        for x in $(ls purgeHaplotigs_04_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done

		### create datander commands
		nblocks=$(getNumOfDbBlocks ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID}/${PROJECT_ID}_CT_M.db)

		### create TANmask commands
        for x in $(seq 1 ${fixblocks})
        do 
            echo "cd ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID} && ${DAZZLER_PATH}/bin/TANmask ${PROJECT_ID}_CT_Z.${x} TAN.${PROJECT_ID}_CT_Z.${x}.las && cd ${myCWD}" 
    	done > mask_${sID}_TANmask_block_${FIX_DB%.db}.${slurmID}.plan
        echo "DAZZLER TANmask $(git --git-dir=${DAZZLER_SOURCE_PATH}/DAMASKER/.git rev-parse --short HEAD)" > mask_${sID}_TANmask_block_${FIX_DB%.db}.${slurmID}.version



        for x in $(seq 1 ${nblocks})
        do 
            echo "cd ${CT_PURGEHAPLOTIGS_OUTDIR}/purgeHaplotigs_${CT_PURGEHAPLOTIGS_RUNID} && PATH=${DAZZLER_PATH}/bin:\${PATH} ${DAZZLER_PATH}/bin/datander${CT_PURGEHAPLOTIGS_DATANDER_OPT} ${PROJECT_ID}_CT_Z.${x} && cd ${myCWD}"
		done > purgeHaplotigs_03_TCdatander_block_${CONT_DB}.${slurmID}.plan
        echo "DAZZLER datander $(git --git-dir=${DAZZLER_SOURCE_PATH}/DAMASKER/.git rev-parse --short HEAD)" > purgeHaplotigs_03_TCdatander_block_${CONT_DB}.${slurmID}.version
 
				
    ### 02_PDminimap2			("01_PDprepInput, 02_PDminimap2, 03_PDcalcuts, 04_PDminimap2, 05_purgedups, 06_statistics")
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