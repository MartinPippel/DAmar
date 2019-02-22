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
	if [[ -z "${RAW_MASH_FASTP_THREADS}" ]]
	then 
		RAW_MASH_FASTP_THREADS=4	
	fi
}

#type-0 [10x - prepare] [1-3]: 01_longrangerBasic, 02_longrangerToScaff10Xinput, 03_bxcheck
#type-1 [10x - de novo] [1-1]: 01_supernova
#type-2 [10x|HiC - kmer-Gsize estimate] [1-2]: 01_jellyfish, 02_genomescope
#type-3 [allData - MASH CONTAMINATION SCREEN] [1-5]: 01_mashPrepare, 02_mashSketch, 03_mashCombine, 04_mashPlot, 05_mashScreen


myTypes=("01_longrangerBasic, 02_longrangerToScaff10Xinput, 03_bxcheck",
"01_supernova", "01_jellyfish, 02_genomescope"
"01_mashPrepare, 02_mashSketch, 03_mashCombine, 04_mashPlot, 05_mashScreen")

#type-0 [10x - prepare] [1-3]: 01_longrangerBasic, 02_longrangerToScaff10Xinput, 03_bxcheck
if [[ ${RAW_QC_TYPE} -eq 0 ]]
then
	### 01_longrangerBasic
    if [[ ${currentStep} -eq 1 ]]
    then
        ### clean up plans 
        for x in $(ls qc_01_*_*_${RAW_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        ## check if 10x data is available
        if [[ -n ${TENX_PATH} && -d "${TENX_PATH}" ]]
        then
        	## we have to trim off the 10x adapter
        	
        	numR1Files=0
			for x in ${TENX_PATH}/${PROJECT_ID}_S*_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz
			do
				if [[ -f ${x} ]]
				then	
					numR1Files=$((${numR1Files}+1))	
				fi
			done
			
			if [[ ${numR1Files} -eq 0 ]]
	        then
	        	(>&2 echo "ERROR - cannot read 10x R1 files with following pattern: ${TENX_PATH}/${PROJECT_ID}_S*_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz")
	        	exit 1
	   		fi
	   		
	   		numR2Files=0
			for x in ${TENX_PATH}/${PROJECT_ID}_S*_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz
			do
				if [[ -f ${x} ]]
				then	
					numR2Files=$((${numR2Files}+1))	
				fi
			done
			
			if [[ ${numR2Files} -eq 0 ]]
	        then
	        	(>&2 echo "ERROR - cannot read 10x R2 files with following pattern: ${TENX_PATH}/${PROJECT_ID}_S*_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz")
	        	exit 1
	   		fi
	   		
	   		if [[ ${numR1Files} -ne ${numR2Files} ]]
	        then
	        	(>&2 echo "ERROR - 10x R1 files ${numR1Files} does not match R2 files ${numR2Files}")
	        	exit 1
	   		fi	   		
	   		
	   		if [[ -d 10x_${PROJECT_ID}_longrangerBasic ]]
	   		then 
	   			echo "mv 10x_${PROJECT_ID}_longrangerBasic 10x_${PROJECT_ID}_longrangerBasic_$(date '+%Y-%m-%d_%H-%M-%S')"
	   		fi
	   		
	   		echo "${LONGRANGER_PATH}/longranger basic --id=10x_${PROJECT_ID}_longrangerBasic --fastqs=${TENX_PATH}"	   					 
		fi > qc_01_longrangerBasic_single_${RAW_DB%.db}.${slurmID}.plan
        
        echo "$(${LONGRANGER_PATH}/longranger basic --version | head -n1)" > qc_01_longrangerBasic_single_${RAW_DB%.db}.${slurmID}.version	
	else
        (>&2 echo "step ${currentStep} in RAW_QC_TYPE ${RAW_QC_TYPE} not supported")
        (>&2 echo "valid steps are: ${myTypes[${RAW_QC_TYPE}]}")
        exit 1            
    fi
#type-1 [10x - de novo] [1-1]: 01_supernova	
elif [[ ${RAW_QC_TYPE} -eq 1 ]]
then 
	### 01_supernova
    if [[ ${currentStep} -eq 1 ]]
    then
        ### clean up plans 
        for x in $(ls qc_01_*_*_${RAW_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        ## check if 10x data is available
        if [[ -n ${TENX_PATH} && -d "${TENX_PATH}" ]]
        then
        	## we have to trim off the 10x adapter
        	
        	numR1Files=0
			for x in ${TENX_PATH}/${PROJECT_ID}_S*_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz
			do
				if [[ -f ${x} ]]
				then	
					numR1Files=$((${numR1Files}+1))	
				fi
			done
			
			if [[ ${numR1Files} -eq 0 ]]
	        then
	        	(>&2 echo "ERROR - cannot read 10x R1 files with following pattern: ${TENX_PATH}/${PROJECT_ID}_S*_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz")
	        	exit 1
	   		fi
	   		
	   		numR2Files=0
			for x in ${TENX_PATH}/${PROJECT_ID}_S*_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz
			do
				if [[ -f ${x} ]]
				then	
					numR2Files=$((${numR2Files}+1))	
				fi
			done
			
			if [[ ${numR2Files} -eq 0 ]]
	        then
	        	(>&2 echo "ERROR - cannot read 10x R2 files with following pattern: ${TENX_PATH}/${PROJECT_ID}_S*_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz")
	        	exit 1
	   		fi
	   		
	   		if [[ ${numR1Files} -ne ${numR2Files} ]]
	        then
	        	(>&2 echo "ERROR - 10x R1 files ${numR1Files} does not match R2 files ${numR2Files}")
	        	exit 1
	   		fi	   		
	   		
	   		if [[ -d 10x_${PROJECT_ID}_supernova ]]
	   		then 
	   			echo "mv 10x_${PROJECT_ID}_supernova 10x_${PROJECT_ID}_supernova_$(date '+%Y-%m-%d_%H-%M-%S')"	   			
	   		fi
	   		
	   		echo "${SUPERNOVA_PATH}/supernova run --id=10x_${PROJECT_ID} --fastqs=${TENX_PATH}"	   					 
		fi > qc_01_supernova_single_${RAW_DB%.db}.${slurmID}.plan
        
        echo "$(${SUPERNOVA_PATH}/supernova run --version | head -n1)" > qc_01_supernova_single_${RAW_DB%.db}.${slurmID}.version
	else
        (>&2 echo "step ${currentStep} in RAW_QC_TYPE ${RAW_QC_TYPE} not supported")
        (>&2 echo "valid steps are: ${myTypes[${RAW_QC_TYPE}]}")
        exit 1            
    fi		
#type-2 [10x|HiC - kmer-Gsize estimate] [1-2]: 01_jellyfish, 02_genomescope
elif [[ ${RAW_QC_TYPE} -eq 2 ]]
then  
	    ### 01_mashPrepare
    if [[ ${currentStep} -eq 1 ]]
    then
        ### clean up plans 
        for x in $(ls qc_01_*_*_${RAW_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
	fi

## mash contamination check
elif [[ ${RAW_QC_TYPE} -eq 3 ]]
then 
    ### 01_mashPrepare
    if [[ ${currentStep} -eq 1 ]]
    then
        ### clean up plans 
        for x in $(ls qc_01_*_*_${RAW_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        doPacbio=0
        do10x=0
        doHic=0
        
        if [[ -n ${DB_PATH} && -d "${DB_PATH}" ]]
        then 
        	for x in ${DB_PATH}/*subreads.bam
        	do
        		if [[ -f $x ]]
        		then 
        			doPacbio=$((${doPacbio}+1))
        			of="pacbio/$(basename ${x%.subreads.bam})_mash.fasta.gz"
        			echo "PATH=${DAZZLER_PATH}/bin:\${PATH} ${DAZZLER_PATH}/bin/dextract -v -f -o $x | gzip > ${of}"
        		fi
        	done
		fi > qc_01_mashPrepare_block_${RAW_DB%.db}.${slurmID}.plan
		
		if [[ -n ${TENX_PATH} && -d "${TENX_PATH}" ]]
        then
        	## we have to trim off the 10x adapter
        	
        	numR1Files=0
			for x in TENX_PATH/${PROJECT_ID}_S*_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz
			do
				if [[ -f ${x} ]]
				then	
					numR1Files=$((${numR1Files}+1))	
				fi
			done
			
			if [[ ${numR1Files} -eq 0 ]]
	        then
	        	(>&2 echo "ERROR - cannot read 10x R1 files with following pattern: ${TENX_PATH}/${PROJECT_ID}_S*_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz")
	        	exit 1
	   		fi
	   		
	   		numR2Files=0
			for x in ${TENX_PATH}/${PROJECT_ID}_S*_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz
			do
				if [[ -f ${x} ]]
				then	
					numR2Files=$((${numR2Files}+1))	
				fi
			done
			
			if [[ ${numR2Files} -eq 0 ]]
	        then
	        	(>&2 echo "ERROR - cannot read 10x R2 files with following pattern: ${TENX_PATH}/${PROJECT_ID}_S*_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz")
	        	exit 1
	   		fi
	   		
	   		if [[ ${numR1Files} -ne ${numR2Files} ]]
	        then
	        	(>&2 echo "ERROR - 10x R1 files ${numR1Files} does not match R2 files ${numR2Files}")
	        	exit 1
	   		fi
	   		
	   		setFastpOptions
	   		
			for r1 in ${TENX_PATH}/${PROJECT_ID}_S[0-9]_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz
			do
				do10x=$((${do10x}+1))
				id=$(dirname ${r1})
				f1=$(basename ${r1})
				f2=$(echo "${f1}" | sed -e "s:_R1_:_R2_:")
				o="${f1%_R1_???.fastq.gz}"											
				
				echo "fastp -i ${id}/${f1} -I ${id}/${f2} -f 23 -G -Q -j 10x/${o}.json -h 10x/${o}.html -w ${RAW_MASH_FASTP_THREADS} -o 10x/${f1} -O 10x/${f2}"				 
			done 
    	fi >> qc_01_mashPrepare_block_${RAW_DB%.db}.${slurmID}.plan
    	
    	if [[ -n ${HIC_PATH} && -d "${HIC_PATH}" ]]
        then
        	numR1Files=0
			for x in ${HIC_PATH}/${PROJECT_ID}_*_*_R1.fastq.gz
			do
				if [[ -f ${x} ]]
				then	
					numR1Files=$((${numR1Files}+1))	
				fi
			done
			
			if [[ ${numR1Files} -eq 0 ]]
	        then
	        	(>&2 echo "ERROR - cannot read HiC R1 files with following pattern: ${HIC_PATH}/${PROJECT_ID}_*_*_R1.fastq.gz")
	        	exit 1
	   		fi
			
			numR2Files=0
			for x in ${HIC_PATH}/${PROJECT_ID}_*_*_R2.fastq.gz
			do
				if [[ -f ${x} ]]
				then	
					numR2Files=$((${numR2Files}+1))	
				fi
			done
			
			if [[ ${numR2Files} -eq 0 ]]
	        then
	        	(>&2 echo "ERROR - cannot read HiC R2 files with following pattern: ${HIC_PATH}/${PROJECT_ID}_*_*_R2.fastq.gz")
	        	exit 1
	   		fi
	   		
	   		if [[ ${numR1Files} -ne ${numR2Files} ]]
	        then
	        	(>&2 echo "ERROR - HiC R1 files ${numR1Files} does not match R2 files ${numR2Files}")
	        	exit 1
	   		fi
	   		
	   		doHic=${numR1Files}
	   		for x in ${SC_HIC_READS}/${PROJECT_ID}_*_*_R[12].fastq.gz
			do
				echo "ln -s -r -f ${x} hic/"  
			done 
    	fi >> qc_01_mashPrepare_block_${RAW_DB%.db}.${slurmID}.plan
    	
    	if [[ ${doPacbio} -gt 0 ]]
    	then 
    		mkdir -p pacbio	
    	fi
    	
    	if [[ ${do10x} -gt 0 ]]
    	then 
    		mkdir -p 10x	
    	fi
    	
    	if [[ ${doHic} -gt 0 ]]
    	then 
    		mkdir -p hic	
    	fi 
    ### 02_mashSketch
    elif [[ ${currentStep} -eq 2 ]]
    then
        ### clean up plans 
        for x in $(ls qc_02_*_*_${RAW_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
    	# check which subdirs are present and contain files
    	# pacbio
    	if [[ -n ${DB_PATH} && -d "${DB_PATH}" && -d pacbio ]]
    	then 
    		for x in pacbio/*_mash.fasta.gz
    		do
    			if [[ -f ${x} ]]
    			then 
    				echo "mash sketch -k 21 -s 10000 -r -m 1 -o ${x%.fasta.gz}.msh ${x}"
    			fi	
    		done	
		fi > qc_02_mashSketch_block_${RAW_DB%.db}.${slurmID}.plan
    
        # 10x
        if [[ -n ${TENX_PATH} && -d "${TENX_PATH}" && -d 10x ]]
    	then 
    		for x in 10x/${PROJECT_ID}_S*_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz
    		do
    			echo "cat ${x} $(echo ${x} | sed -e "s:_R1_:_R2_:") mash sketch -k 21 -s 10000 -r -m 2 -o $(echo ${x%.fastq.gz}.msh | sed -e "s:_R1::") -"
    		done
    	fi >> qc_02_mashSketch_block_${RAW_DB%.db}.${slurmID}.plan
    	
    	# hic
    	if [[ -n ${HIC_PATH} && -d "${HIC_PATH}" && -d hic ]]
    	then
    		for x in ${HIC_PATH}/${PROJECT_ID}_*_*_R1.fastq.gz
			do
				echo "cat ${x} ${x%_R1.fastq.gz}_R2.fastq.gz | mash sketch -k 21 -s 10000 -r -m 2 -o ${x%_R1.fastq.gz}.msh -"
			done    		
    	fi >> qc_02_mashSketch_block_${RAW_DB%.db}.${slurmID}.plan
    	
    	echo "mash $(${PACBIO_BASE_ENV} && mash --version && ${PACBIO_BASE_ENV_DEACT})" > qc_02_mashSketch_block_${RAW_DB%.db}.${slurmID}.version
    ### 03_mashCombine
    elif [[ ${currentStep} -eq 3 ]]
    then
    	### clean up plans 
        for x in $(ls qc_03_*_*_${RAW_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        if [[ -d pacbio ]]
        then 
        	ls pacbio/*.msh > ${PROJECT_ID}_mash.files
    	fi
    	
    	if [[ -d 10x ]]
        then 
        	ls 10x/*.msh >> ${PROJECT_ID}_mash.files
    	fi 
    	
    	if [[ -d hic ]]
        then 
        	ls hic/*.msh >> ${PROJECT_ID}_mash.files
    	fi
    	        
        echo "mash paste -l ${PROJECT_ID}.msh ${PROJECT_ID}_mash.files" > qc_03_mashCombine_single_${RAW_DB%.db}.${slurmID}.plan
        echo "mash dist -t ${PROJECT_ID}.msh ${PROJECT_ID}.msh > ${PROJECT_ID}.tbl" >> qc_03_mashCombine_single_${RAW_DB%.db}.${slurmID}.plan
        echo "mash $(${PACBIO_BASE_ENV} && mash --version && ${PACBIO_BASE_ENV_DEACT})" > qc_03_mashCombine_single_${RAW_DB%.db}.${slurmID}.version
	### 04_mashPlot
    elif [[ ${currentStep} -eq 4 ]]
    then
    	### clean up plans 
        for x in $(ls qc_04_*_*_${RAW_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
    
    	echo "head -n 1 ${PROJECT_ID}.tbl |awk '{for (i=2; i <=NF; i++) print \$i}' |awk -F "/" '{print \$NF}' |sed s/.subreads.fast[aq].gz//g |sed s/.fast[aq].gz//g |sed s/.fast[aq]//g > ${PROJECT_ID}.key" > qc_04_mashPlot_single_${RAW_DB%.db}.${slurmID}.plan
    	echo "Rscript ${MARVEL_PATH}/scripts/mashPlot.R ${PROJECT_ID}" >> qc_04_mashPlot_single_${RAW_DB%.db}.${slurmID}.plan
    	echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > qc_04_mashPlot_single_${RAW_DB%.db}.${slurmID}.version
    ### 05_mashScreen
    elif [[ ${currentStep} -eq 5 ]]
    then
    	### clean up plans 
        for x in $(ls qc_05_*_*_${RAW_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        if [[ -z ${MASH_REF_GENOMES} || ! -f ${MASH_REF_GENOMES} ]]
        then
        	(>&2 echo "[ERROR] mash pipeline step 05 mashScreen failed. Set variable MASH_REF_GENOMES to an appropriate reference genome msh file")
        	exit 1	
    	fi
        	
        mkdir -p screen 
        
        threads=24
		if [[ "${SLURM_PARTITION}" == "gpu" ]]
		then 
			threads=40
		elif [[ "${SLURM_PARTITION}" == "bigmem" ]]
		then 
			threads=48
		fi
        
        for x in pacbio/*.fasta.gz 10x/*.fastq.gz hic/*.fastq.gz
        do
        	if [[ -f ${x} ]]
        	then
        		out=screen/$(basename ${x%.fast[aq].gz}).conta
        		echo "mash screen -p ${threads} -w ${MASH_REF_GENOMES} ${x} > ${out}"
        	fi
		done > qc_05_mashScreen_block_${RAW_DB%.db}.${slurmID}.plan   
		echo "mash $(${PACBIO_BASE_ENV} && mash --version && ${PACBIO_BASE_ENV_DEACT})" > qc_05_mashScreen_block_${RAW_DB%.db}.${slurmID}.version     
    else
        (>&2 echo "step ${currentStep} in RAW_QC_TYPE ${RAW_QC_TYPE} not supported")
        (>&2 echo "valid steps are: ${myTypes[${RAW_QC_TYPE}]}")
        exit 1            
    fi        
else
    (>&2 echo "unknown RAW_QC_TYPE ${RAW_QC_TYPE}")
    (>&2 echo "supported types")
    x=0; while [ $x -lt ${#myTypes[*]} ]; do (>&2 echo "${myTypes[${x}]}"); done 
    exit 1
fi

exit 0