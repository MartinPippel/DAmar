#!/bin/bash 

configFile=$1
currentStep=$2
slurmID=$3

if [[ ! -f ${configFile} ]]
then 
    (>&2 echo "cannot access config file ${configFile}")
    exit 1
fi

source ${configFile}

if [[ ! -n "${RAW_DASCOVER_TYPE}" ]]
then 
    (>&2 echo "cannot create repmask jobs if varibale RAW_REPMASK_TYPE is not set.")
    exit 1
fi

if [[ ! -n ${RAW_DAZZ_DB} ]]
then 
    (>&2 echo "raw database unknown - You have to set the variable RAW_DAZZLER_DB")
    exit 1
fi

if [[ ! -f ${RAW_DAZZ_DB%.db}.db ]]
then 
    (>&2 echo "raw database ${RAW_DAZZ_DB%.db}.db missing")
    exit 1
fi

function getNumOfDbBlocks()
{
    if [[ ! -f ${RAW_DAZZ_DB%.db}.db ]]
    then
        (>&2 echo "raw database ${RAW_DAZZ_DB%.db}.db not found")
        exit 1
    fi

    blocks=$(grep block ${RAW_DAZZ_DB%.db}.db | awk '{print $3}')
    if [[ ! -n $blocks ]]
    then 
        (>&2 echo "raw database ${RAW_DAZZ_DB%.db}.db has not been partitioned. Run DBsplit first!")
        exit 1
    fi 
    echo ${blocks}
}

nblocks=$(getNumOfDbBlocks)

function setDBdustOptions()
{
    DASCOVER_DBDUST_OPT=""
    if [[ -n ${RAW_DASCOVER_DBDUST_BIAS} && ${RAW_DASCOVER_DBDUST_BIAS} -ge 1 ]]
    then
        DASCOVER_DBDUST_OPT="${DASCOVER_DBDUST_OPT} -b"
    fi
}

function setDatanderOptions()
{
    ### find and set datander options 
    DASCOVER_DATANDER_OPT=""
    if [[ -n ${RAW_DASCOVER_DATANDER_THREADS} ]]
    then
        DASCOVER_DATANDER_OPT="${DASCOVER_DATANDER_OPT} -T${RAW_DASCOVER_DATANDER_THREADS}"
    fi
    if [[ -n ${RAW_DASCOVER_DATANDER_MINLEN} ]]
    then
        DASCOVER_DATANDER_OPT="${DASCOVER_DATANDER_OPT} -l${RAW_DASCOVER_DATANDER_MINLEN}"
    fi
}

function setLAcheckOptions()
{
	# Check that .las is in sorted order BY DEFAULT
	DASCOVER_LACHECK_OPT=" -S"
	
	if [[ -z ${DASCOVER_LACHECK_BLOCKCMP} ]]
	then 
		DASCOVER_LACHECK_BLOCKCMP=4	
	fi
	if [[ -n ${DASCOVER_LACHECK_VERBOSE} ]]
	then 
		DASCOVER_LACHECK_OPT="${DASCOVER_LACHECK_OPT} -v"
	fi  
}

function setTANmaskOptions()
{
	DASCOVER_TANMASK_OPT=""
	
	if [[ -z ${DASCOVER_TANMASK_BLOCKCMP} ]]
	then 
		DASCOVER_TANMASK_BLOCKCMP=4	
	fi
	if [[ -n ${DASCOVER_TANMASK_VERBOSE} ]]
	then 
		DASCOVER_TANMASK_OPT="${DASCOVER_TANMASK_OPT} -v"
	fi  
}

function setDaligerOptions()
{
    DASCOVER_DALIGNER_OPT=""
    if [[ -n ${RAW_DASCOVER_DALIGNER_KMER} && ${RAW_DASCOVER_DALIGNER_KMER} -gt 0 ]]
    then
        DASCOVER_DALIGNER_OPT="${DASCOVER_DALIGNER_OPT} -k${RAW_DASCOVER_DALIGNER_KMER}"
    fi
    if [[ -n ${RAW_DASCOVER_DALIGNER_ERR} ]]
    then
        DASCOVER_DALIGNER_OPT="${DASCOVER_DALIGNER_OPT} -e${RAW_DASCOVER_DALIGNER_ERR}"
    fi
    if [[ -n ${RAW_DASCOVER_DALIGNER_BIAS} && ${RAW_DASCOVER_DALIGNER_BIAS} -eq 1 ]]
    then
        DASCOVER_DALIGNER_OPT="${DASCOVER_DALIGNER_OPT} -b"
    fi
    if [[ -n ${RAW_DASCOVER_DALIGNER_OLEN} ]]
    then
        DASCOVER_DALIGNER_OPT="${DASCOVER_DALIGNER_OPT} -l${RAW_DASCOVER_DALIGNER_OLEN}"
    fi    
    if [[ -n ${RAW_DASCOVER_DALIGNER_MEM} && ${RAW_DASCOVER_DALIGNER_MEM} -gt 0 ]]
    then
        DASCOVER_DALIGNER_OPT="${DASCOVER_DALIGNER_OPT} -M${RAW_DASCOVER_DALIGNER_MEM}"
    fi    
    if [[ -n ${RAW_DASCOVER_DALIGNER_HITS} ]]
    then
        DASCOVER_DALIGNER_OPT="${DASCOVER_DALIGNER_OPT} -h${RAW_DASCOVER_DALIGNER_HITS}"
    fi        
    if [[ -n ${RAW_DASCOVER_DALIGNER_T} ]]
    then
        DASCOVER_DALIGNER_OPT="${DASCOVER_DALIGNER_OPT} -t${RAW_DASCOVER_DALIGNER_T}"
    fi  
    if [[ -n ${RAW_DASCOVER_DALIGNER_MASK} ]]
    then
        for x in ${RAW_DASCOVER_DALIGNER_MASK}
        do 
            DASCOVER_DALIGNER_OPT="${DASCOVER_DALIGNER_OPT} -m ${x}"
        done
    fi
    if [[ -n ${RAW_DASCOVER_DALIGNER_TRACESPACE} && ${RAW_DASCOVER_DALIGNER_TRACESPACE} -gt 0 ]]
    then
        DASCOVER_DALIGNER_OPT="${DASCOVER_DALIGNER_OPT} -s${RAW_DASCOVER_DALIGNER_TRACESPACE}"
    fi
    
    if [[ -n ${THREADS_daligner} ]]
    then 
        DASCOVER_DALIGNER_OPT="${DASCOVER_DALIGNER_OPT} -T${THREADS_daligner}"
    fi
    
    if [[ -z ${DASCOVER_DALIGNER_BLOCKCMP} ]]
    then
    	DASCOVER_DALIGNER_BLOCKCMP=4	
	fi
	
	if [[ -z ${RAW_DASCOVER_DALIGNER_FORBLOCK} ]]
	then
		RAW_DASCOVER_DALIGNER_FORBLOCK=1	
	fi
}

function setREPmaskOptions()
{
	DASCOVER_REPMASK_OPT=""
	
	if [[ -z ${DASCOVER_REPMASK_BLOCKCMP} ]]
	then 
		DASCOVER_REPMASK_BLOCKCMP=4	
	fi
	
	if [[ -n ${DASCOVER_REPMASK_VERBOSE} ]]
	then 
		DASCOVER_REPMASK_OPT="${DASCOVER_REPMASK_OPT} -v"
	fi
	
	if [[ -n ${DASCOVER_REPMASK_COVERAGE} && ${DASCOVER_REPMASK_COVERAGE} -gt 0  ]]
	then 
		DASCOVER_REPMASK_OPT="${DASCOVER_REPMASK_OPT} -c${DASCOVER_REPMASK_COVERAGE}"
	else
		if [[ -z ${GSIZE} ]]
		then
			(>&2 echo "REPmask needs a coverage estimation! You have to set either DASCOVER_REPMASK_COVERAGE or a rough estimate of the genome size: GSIZE")
        	exit 1	
		fi
		
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
		
		blockSize=$(grep size ${RAW_DAZZ_DB%.db}.db  | awk '{print $3}')
		
		if [[ ${blockSize} -le ${gsize} ]]
		then 
			DASCOVER_REPMASK_OPT="${DASCOVER_REPMASK_OPT} -c10"
		else
			cov=$(echo "scale=0;((${blockSize}/${gsize})+1)*10" |bc)
			DASCOVER_REPMASK_OPT="${DASCOVER_REPMASK_OPT} -c${cov}"
		fi
	fi	
}

function setLAmergeOptions()
{
    DASCOVER_LAMERGE_OPT=""
    if [[ -n ${RAW_DASCOVER_LAMERGE_VERBOSE} && ${RAW_DASCOVER_LAMERGE_VERBOSE} -ge 1 ]]
    then
        DASCOVER_LAMERGE_OPT="${DASCOVER_LAMERGE_OPT} -v"
    fi
}

function setDAScoverOptions()
{		
	DASCOVER_DASCOVER_OPT=""
    if [[ -n ${RAW_DASCOVER_DASCOVER_VERBOSE} && ${RAW_DASCOVER_DASCOVER_VERBOSE} -ge 1 ]]
    then
        DASCOVER_DASCOVER_OPT="${DASCOVER_DASCOVER_OPT} -v"
    fi
    
    if [[ -n ${RAW_DASCOVER_DASCOVER_REPEAT} ]]
    then
    	for x in ${RAW_DASCOVER_DASCOVER_REPEAT}
    	do
    		DASCOVER_DASCOVER_OPT="${DASCOVER_DASCOVER_OPT} -m${x}"	
    	done        
    fi
    
	if [[ -n ${RAW_DASCOVER_DASCOVER_HGAPLEN} && ${RAW_DASCOVER_DASCOVER_HGAPLEN} -ge 1 ]]
    then
        DASCOVER_DASCOVER_OPT="${DASCOVER_DASCOVER_OPT} -H${RAW_DASCOVER_DASCOVER_HGAPLEN}"
    fi
}

## ensure some paths

if [[ -z "${DAZZLER_SOURCE_PATH}" ]]
then 
    (>&2 echo "ERROR - You have to set DAZZLER_SOURCE_PATH. Used to report git version.")
    exit 1
fi

myTypes=("1-DBdust, 2-datander, 3-TANmask, 4-CheckTan, 5-rmTan, 6-daligner, 7-CheckDaligner, 8-REPmask, 9-rmDaligner, 10-daligner, 11-CheckDaligner, 12-LAmerge, 13-CheckMerge, 14-rmDaligner, 15-DAScover, 16-REPcover")
# type-0 steps [1-13]: 1-DBdust, 2-datander, 3-CheckTan, 4-TANmask, 5-rmTan, 6-daligner, 7-CheckDaligner, 8-REPmask, 9-rmDaligner, 10-daligner, 11-CheckDaligner, 12-LAmerge, 13-CheckMerge, 14-rmDaligner, 15-DAScover, 16-REPcover
if [[ ${RAW_DASCOVER_TYPE} -eq 0 ]]
then 
    if [[ ${currentStep} -eq 1 ]]
    then
        ### clean up plans 
        for x in $(ls cover_01_*_*_${RAW_DAZZ_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        ### find and set DBdust options 
        setDBdustOptions
        ### create DBdust commands 
        for x in $(seq 1 ${nblocks})
        do 
            echo "${DAZZLER_PATH}/bin/DBdust${DASCOVER_DBDUST_OPT} ${RAW_DAZZ_DB%.db}.${x}"
		done > cover_01_DBdust_block_${RAW_DAZZ_DB%.db}.${slurmID}.plan
        echo "DAZZLER $(git --git-dir=${DAZZLER_SOURCE_PATH}/DAZZ_DB/.git rev-parse --short HEAD)" > cover_01_DBdust_block_${RAW_DAZZ_DB%.db}.${slurmID}.version
    elif [[ ${currentStep} -eq 2 ]]
    then 
        ### clean up plans 
        for x in $(ls cover_02_*_*_${RAW_DAZZ_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done     
        ### find and set datander options 
        setDatanderOptions
        ### create datander commands
        for x in $(seq 1 ${nblocks})
        do 
            echo "PATH=${DAZZLER_PATH}/bin:\${PATH} ${DAZZLER_PATH}/bin/datander${DASCOVER_DATANDER_OPT} ${RAW_DAZZ_DB%.db}.${x}"
		done > cover_02_datander_block_${RAW_DAZZ_DB%.db}.${slurmID}.plan
		echo "DAMASKER $(git --git-dir=${DAZZLER_SOURCE_PATH}/DAMASKER/.git rev-parse --short HEAD)" > cover_02_datander_block_${RAW_DAZZ_DB%.db}.${slurmID}.version        
    elif [[ ${currentStep} -eq 3 ]]
    then 
        ### clean up plans 
        for x in $(ls cover_03_*_*_${RAW_DAZZ_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done     
        ### find and set LAcheck options 
        setLAcheckOptions
        ### create LAcheck commands
        x=1;
        while [[ ${x} -le ${nblocks} ]]
        do 
        	y=$((${x}+${DASCOVER_LACHECK_BLOCKCMP}))
        	if [[ $y -gt ${nblocks} ]]
        	then
        		y=${nblocks}	
        	fi 
        	echo "${DAZZLER_PATH}/bin/LAcheck${DASCOVER_LACHECK_OPT} ${RAW_DAZZ_DB%.db} TAN.${RAW_DAZZ_DB%.db}.@${x}-${y} || exit 1"
        	x=$((${y}+1)) 
		done > cover_03_LAcheck_block_${RAW_DAZZ_DB%.db}.${slurmID}.plan
        echo "DALIGNER $(git --git-dir=${DAZZLER_SOURCE_PATH}/DALIGNER/.git rev-parse --short HEAD)" > cover_03_LAcheck_block_${RAW_DAZZ_DB%.db}.${slurmID}.version
    elif [[ ${currentStep} -eq 4 ]]
    then 
        ### clean up plans 
        for x in $(ls cover_04_*_*_${RAW_DAZZ_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done     
        ### find and set TANmask options 
        setTANmaskOptions
        ### create TANmask commands
        x=1;
        while [[ ${x} -le ${nblocks} ]]
        do 
        	y=$((${x}+${DASCOVER_TANMASK_BLOCKCMP}))
        	if [[ $y -gt ${nblocks} ]]
        	then
        		y=${nblocks}	
        	fi
        	echo "${DAZZLER_PATH}/bin/TANmask${DASCOVER_TANMASK_OPT} ${RAW_DAZZ_DB%.db} TAN.${RAW_DAZZ_DB%.db}.@${x}-${y}"
        	x=$((${y}+1)) 
		done > cover_04_TANmask_block_${RAW_DAZZ_DB%.db}.${slurmID}.plan
        echo "DAMASKER $(git --git-dir=${DAZZLER_SOURCE_PATH}/DAMASKER/.git rev-parse --short HEAD)" > cover_04_TANmask_block_${RAW_DAZZ_DB%.db}.${slurmID}.version    
    elif [[ ${currentStep} -eq 5 ]]
    then 
        ### clean up plans 
        for x in $(ls cover_05_*_*_${RAW_DAZZ_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        ### create rmTAN command
        x=1;
        while [[ ${x} -le ${nblocks} ]]
        do
        	flist=""
        	NUM_FILES_TO_DEL=50
        	y=$((${x}+${NUM_FILES_TO_DEL}))
        	if [[ ${y} -gt ${nblocks} ]]
        	then
        		y=${nblocks}
        	fi
        	for z in $(seq ${x} ${y})
        	do
        		flist="${flist} TAN.${RAW_DAZZ_DB%.db}.${z}.las"
        	done
        	x=$((${y}+1))
        	echo "rm${flist}"
		done  > cover_05_rmTan_single_${RAW_DAZZ_DB%.db}.${slurmID}.plan            
    elif [[ ${currentStep} -eq 6 ]]
    then
        for x in $(ls cover_06_*_*_${RAW_DAZZ_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        ### find and set daligner options 
        setDaligerOptions

        ### create daligner commands
        for x in $(seq 1 ${nblocks})
        do
    		if [[ -n ${RAW_DASCOVER_DALIGNER_NUMACTL} && ${RAW_DASCOVER_DALIGNER_NUMACTL} -gt 0 ]] && [[ "x${SLURM_NUMACTL}" == "x" || ${SLURM_NUMACTL} -eq 0 ]]
            then
                if [[ $((${x} % 2)) -eq  0 ]]
                then
                    NUMACTL="numactl -m0 -N0 "
                else
                    NUMACTL="numactl -m1 -N1 "    
                fi
            else
                NUMACTL=""
            fi
            echo "PATH=${DAZZLER_PATH}/bin:\${PATH} ${NUMACTL}${DAZZLER_PATH}/bin/daligner${DASCOVER_DALIGNER_OPT} -mdust -mtan ${RAW_DAZZ_DB%.db}.${x} ${RAW_DAZZ_DB%.db}.${x}"
		done > cover_06_daligner_block_${RAW_DAZZ_DB%.db}.${slurmID}.plan
        echo "DALIGNER $(git --git-dir=${DAZZLER_SOURCE_PATH}/DALIGNER/.git rev-parse --short HEAD)" > cover_06_daligner_block_${RAW_DAZZ_DB%.db}.${slurmID}.version
	elif [[ ${currentStep} -eq 7 ]]
    then 
        ### clean up plans 
        for x in $(ls cover_07_*_*_${RAW_DAZZ_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done     
        ### find and set LAcheck options 
        setLAcheckOptions
        ### create LAcheck commands
        x=1;
        while [[ ${x} -le ${nblocks} ]]
        do
        	flist=""
        	NUM_FILES_TO_DEL=${DASCOVER_LACHECK_BLOCKCMP}
        	y=$((${x}+${NUM_FILES_TO_DEL}))
        	if [[ ${y} -gt ${nblocks} ]]
        	then
        		y=${nblocks}
        	fi
        	for z in $(seq ${x} ${y})
        	do
        		flist="${flist} ${RAW_DAZZ_DB%.db}.${z}.${RAW_DAZZ_DB%.db}.${z}"
        	done
        	echo "${DAZZLER_PATH}/bin/LAcheck${DASCOVER_LACHECK_OPT} ${RAW_DAZZ_DB%.db}${flist} || exit 1"
        	x=$((${y}+1))
		done > cover_07_LAcheck_block_${RAW_DAZZ_DB%.db}.${slurmID}.plan
        echo "DALIGNER $(git --git-dir=${DAZZLER_SOURCE_PATH}/DALIGNER/.git rev-parse --short HEAD)" > cover_07_LAcheck_block_${RAW_DAZZ_DB%.db}.${slurmID}.version
	elif [[ ${currentStep} -eq 8 ]]
    then 
        ### clean up plans 
        for x in $(ls cover_08_*_*_${RAW_DAZZ_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done     
        ### find and set LAcheck options 
        setREPmaskOptions
        ### create REPmask commands
        x=1;
        while [[ ${x} -le ${nblocks} ]]
        do
        	flist=""
        	NUM_FILES_TO_DEL=${DASCOVER_REPMASK_BLOCKCMP}
        	y=$((${x}+${NUM_FILES_TO_DEL}))
        	if [[ ${y} -gt ${nblocks} ]]
        	then
        		y=${nblocks}
        	fi
        	for z in $(seq ${x} ${y})
        	do
        		flist="${flist} ${RAW_DAZZ_DB%.db}.${z}.${RAW_DAZZ_DB%.db}.${z}"
        	done
        	echo "${DAZZLER_PATH}/bin/REPmask${DASCOVER_REPMASK_OPT} ${RAW_DAZZ_DB%.db}${flist} || exit 1"
        	x=$((${y}+1))
		done > cover_08_REPmask_block_${RAW_DAZZ_DB%.db}.${slurmID}.plan
        echo "DAMASKER $(git --git-dir=${DAZZLER_SOURCE_PATH}/DAMASKER/.git rev-parse --short HEAD)" > cover_08_REPmask_block_${RAW_DAZZ_DB%.db}.${slurmID}.version
	elif [[ ${currentStep} -eq 9 ]]
    then 
        ### clean up plans 
        for x in $(ls cover_09_*_*_${RAW_DAZZ_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        ### create rmDaligner command
        x=1;
        while [[ ${x} -le ${nblocks} ]]
        do
        	flist=""
        	NUM_FILES_TO_DEL=50
        	y=$((${x}+${NUM_FILES_TO_DEL}))
        	if [[ ${y} -gt ${nblocks} ]]
        	then
        		y=${nblocks}
        	fi
        	for z in $(seq ${x} ${y})
        	do
        		flist="${flist} ${RAW_DAZZ_DB%.db}.${z}.${RAW_DAZZ_DB%.db}.${z}.las"
        	done
        	x=$((${y}+1))
        	echo "rm${flist}"
		done  > cover_09_rmDaligner_single_${RAW_DAZZ_DB%.db}.${slurmID}.plan  
	elif [[ ${currentStep} -eq 10 ]]
    then
        for x in $(ls cover_10_*_*_${RAW_DAZZ_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        ### find and set daligner options 
        setDaligerOptions

        ### create daligner commands
       	if [[ -n ${RAW_DASCOVER_DALIGNER_NUMACTL} && ${RAW_DASCOVER_DALIGNER_NUMACTL} -gt 0 ]] && [[ "x${SLURM_NUMACTL}" == "x" || ${SLURM_NUMACTL} -eq 0 ]]
		then
        	NUMACTL="numactl -m1 -N1 "    
       	else
        	NUMACTL=""
        fi
        echo "PATH=${DAZZLER_PATH}/bin:\${PATH} ${NUMACTL}${DAZZLER_PATH}/bin/daligner${DASCOVER_DALIGNER_OPT} -mdust -mtan -mrep ${RAW_DAZZ_DB%.db}.${RAW_DASCOVER_DALIGNER_FORBLOCK} ${RAW_DAZZ_DB%.db}.${RAW_DASCOVER_DALIGNER_FORBLOCK}" > cover_10_daligner_block_${RAW_DAZZ_DB%.db}.${slurmID}.plan
        x=1;
        while [[ ${x} -le ${nblocks} ]]
        do
            if [[ -n ${RAW_DASCOVER_DALIGNER_NUMACTL} && ${RAW_DASCOVER_DALIGNER_NUMACTL} -gt 0 ]] && [[ "x${SLURM_NUMACTL}" == "x" || ${SLURM_NUMACTL} -eq 0 ]]
            then
                if [[ $((${x} % 2)) -eq  0 ]]
                then
                    NUMACTL="numactl -m0 -N0 "
                else
                    NUMACTL="numactl -m1 -N1 "    
                fi
            else
                NUMACTL=""
            fi
            flist=""
            y=$((${x}+${DASCOVER_DALIGNER_BLOCKCMP}))
            if [[ $y -gt ${nblocks} ]]
            then
            	y=${nblocks}
        	fi
            for z in $(seq ${x} ${y})
            do
            	if [[ ${z} -eq ${RAW_DASCOVER_DALIGNER_FORBLOCK} ]]
            	then
            		continue
            	fi 
            	flist="${flist} ${RAW_DAZZ_DB%.db}.${z}"
            done
            if [[ -n "${flist}" ]]
            then
            	echo "PATH=${DAZZLER_PATH}/bin:\${PATH} ${NUMACTL}${DAZZLER_PATH}/bin/daligner${DASCOVER_DALIGNER_OPT} -mdust -mtan -mrep -A ${RAW_DAZZ_DB%.db}.${RAW_DASCOVER_DALIGNER_FORBLOCK} ${flist}"
        	fi
        	x=$((${y}+1))
		done >> cover_10_daligner_block_${RAW_DAZZ_DB%.db}.${slurmID}.plan
        echo "DALIGNER $(git --git-dir=${DAZZLER_SOURCE_PATH}/DALIGNER/.git rev-parse --short HEAD)" > cover_10_daligner_block_${RAW_DAZZ_DB%.db}.${slurmID}.version
	elif [[ ${currentStep} -eq 11 ]]
    then 
        ### clean up plans 
        for x in $(ls cover_11_*_*_${RAW_DAZZ_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done     
        ### find and set LAcheck options 
        setLAcheckOptions
        ### create LAcheck commands
        x=1;
        while [[ ${x} -le ${nblocks} ]]
        do 
        	y=$((${x}+${DASCOVER_LACHECK_BLOCKCMP}))
        	if [[ $y -gt ${nblocks} ]]
        	then
        		y=${nblocks}	
        	fi 
        	echo "${DAZZLER_PATH}/bin/LAcheck${DASCOVER_LACHECK_OPT} ${RAW_DAZZ_DB%.db} ${RAW_DAZZ_DB%.db}.${RAW_DASCOVER_DALIGNER_FORBLOCK}.${RAW_DAZZ_DB%.db}.@${x}-${y} || exit 1"
        	x=$((${y}+1)) 
		done > cover_11_LAcheck_block_${RAW_DAZZ_DB%.db}.${slurmID}.plan
        echo "DALIGNER $(git --git-dir=${DAZZLER_SOURCE_PATH}/DALIGNER/.git rev-parse --short HEAD)" > cover_11_LAcheck_block_${RAW_DAZZ_DB%.db}.${slurmID}.version   
	elif [[ ${currentStep} -eq 12 ]]
    then 
        ### clean up plans 
        for x in $(ls cover_12_*_*_${RAW_DAZZ_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done     
        ### find and set LAcheck options 
        setLAmergeOptions
        ### create LAmerge command
		echo "PATH=${DAZZLER_PATH}/bin:\${PATH} ${DAZZLER_PATH}/bin/LAmerge${DASCOVER_LAMERGE_OPT} ${RAW_DAZZ_DB%.db}.${nblocks}.${RAW_DASCOVER_DALIGNER_FORBLOCK} ${RAW_DAZZ_DB%.db}.${RAW_DASCOVER_DALIGNER_FORBLOCK}.${RAW_DAZZ_DB%.db}.@1-${nblocks} || exit 1" > cover_12_LAmerge_single_${RAW_DAZZ_DB%.db}.${slurmID}.plan
        echo "DALIGNER $(git --git-dir=${DAZZLER_SOURCE_PATH}/DALIGNER/.git rev-parse --short HEAD)" > cover_12_LAmerge_single_${RAW_DAZZ_DB%.db}.${slurmID}.version
	elif [[ ${currentStep} -eq 13 ]]
    then 
        ### clean up plans 
        for x in $(ls cover_13_*_*_${RAW_DAZZ_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done     
        ### find and set LAcheck options 
        setLAcheckOptions
        ### create LAmerge command
		echo "${DAZZLER_PATH}/bin/LAcheck${DASCOVER_LACHECK_OPT} ${RAW_DAZZ_DB%.db} ${RAW_DAZZ_DB%.db}.${nblocks}.${RAW_DASCOVER_DALIGNER_FORBLOCK} || exit 1" > cover_13_LAcheck_single_${RAW_DAZZ_DB%.db}.${slurmID}.plan
        echo "DALIGNER $(git --git-dir=${DAZZLER_SOURCE_PATH}/DALIGNER/.git rev-parse --short HEAD)" > cover_13_LAcheck_single_${RAW_DAZZ_DB%.db}.${slurmID}.version
	elif [[ ${currentStep} -eq 14 ]]
    then 
        ### clean up plans 
        for x in $(ls cover_14_*_*_${RAW_DAZZ_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        ### create rmDaligner command
        x=1;
        while [[ ${x} -le ${nblocks} ]]
        do
        	flist=""
        	NUM_FILES_TO_DEL=50
        	y=$((${x}+${NUM_FILES_TO_DEL}))
        	if [[ ${y} -gt ${nblocks} ]]
        	then
        		y=${nblocks}
        	fi
        	for z in $(seq ${x} ${y})
        	do
        		flist="${flist} ${RAW_DAZZ_DB%.db}.${RAW_DASCOVER_DALIGNER_FORBLOCK}.${RAW_DAZZ_DB%.db}.${z}.las"
        	done
        	echo "rm${flist}"
        	x=$((${y}+1))
		done  > cover_14_rmDaligner_single_${RAW_DAZZ_DB%.db}.${slurmID}.plan
	elif [[ ${currentStep} -eq 15 ]]
    then 
        ### clean up plans 
        for x in $(ls cover_15_*_*_${RAW_DAZZ_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done     
        ### find and set LAcheck options 
        setDAScoverOptions
        ### create DAScover command
		echo "${DAZZLER_PATH}/bin/DAScover${DASCOVER_DASCOVER_OPT} ${RAW_DAZZ_DB%.db} ${RAW_DAZZ_DB%.db}.${nblocks}.${RAW_DASCOVER_DALIGNER_FORBLOCK}" > cover_15_DAScover_single_${RAW_DAZZ_DB%.db}.${slurmID}.plan
        echo "DASCRUBBER $(git --git-dir=${DAZZLER_SOURCE_PATH}/DASCRUBBER/.git rev-parse --short HEAD)" > cover_15_DAScover_single_${RAW_DAZZ_DB%.db}.${slurmID}.version
	elif [[ ${currentStep} -eq 16 ]]
    then
        ### clean up plans 
        for x in $(ls mask_16_*_*_${RAW_DAZZ_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
	
        ### create REPcover commands 
		echo "${DAZZLER_PATH}/bin/REPcover ${RAW_DAZZ_DB%.db}.${RAW_DASCOVER_DALIGNER_FORBLOCK} | tee effectiveCov_${RAW_DAZZ_DB%.db}_${slurmID}_forBlock${RAW_DASCOVER_DALIGNER_FORBLOCK}.txt" > cover_16_REPcover_single_${RAW_DAZZ_DB%.db}.${slurmID}.plan
        echo "DASCRUBBER $(git --git-dir=${DAZZLER_SOURCE_PATH}/DASCRUBBER/.git rev-parse --short HEAD)" > cover_16_DAScover_single_${RAW_DAZZ_DB%.db}.${slurmID}.version
    else 
        (>&2 echo "step ${currentStep} in RAW_DASCOVER_TYPE ${RAW_DASCOVER_TYPE} not supported")
        (>&2 echo "valid steps are: ${myTypes[${RAW_DASCOVER_TYPE}]}")
        exit 1        
    fi    
else
    (>&2 echo "unknown RAW_DASCOVER_TYPE ${RAW_DASCOVER_TYPE}")
    (>&2 echo "supported types")
    x=0; while [ $x -lt ${#myTypes[*]} ]; do (>&2 echo "${myTypes[${x}]}"); done
    exit 1
fi

exit 0