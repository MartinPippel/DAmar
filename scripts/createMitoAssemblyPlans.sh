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

function getSubDirName()
{
    runID=$1
    blockID=$2

    dname="d${runID}"

    if [[ $runID -lt 10 ]]
    then 
        dname="d00${runID}"
    elif [[ $runID -lt 100 ]]
    then 
        dname="d0${runID}"
    fi

    bname="${blockID}"

    if [[ ${blockID} -lt 10 ]]
    then 
        bname="0000${blockID}"
    elif [[ ${blockID} -lt 100 ]]
    then 
        bname="000${blockID}"
    elif [[ ${blockID} -lt 1000 ]]
    then 
        bname="00${blockID}"           
    elif [[ ${blockID} -lt 10000 ]]
    then 
        bname="0${blockID}"           
    fi
    echo ${dname}_${bname}                 
}

rawblocks=0

if [[ -f ${RAW_DB%.db}.db ]]
then 
	rawblocks=$(getNumOfDbBlocks ${RAW_DB%.db}.db)	
fi

function setDaligerOptions()
{
    MITO_DALIGNER_OPT=""
    if [[ -n ${RAW_MITO_DALIGNER_IDENTITY_OVLS} && ${RAW_MITO_DALIGNER_IDENTITY_OVLS} -gt 0 ]]
    then
        MITO_DALIGNER_OPT="${MITO_DALIGNER_OPT} -I"
    fi
    if [[ -n ${RAW_MITO_DALIGNER_KMER} && ${RAW_MITO_DALIGNER_KMER} -gt 0 ]]
    then
        MITO_DALIGNER_OPT="${MITO_DALIGNER_OPT} -k ${RAW_MITO_DALIGNER_KMER}"
    fi
    if [[ -n ${RAW_MITO_DALIGNER_ERR} ]]
    then
        MITO_DALIGNER_OPT="${MITO_DALIGNER_OPT} -e ${RAW_MITO_DALIGNER_ERR}"
    fi
    if [[ -n ${RAW_MITO_DALIGNER_BIAS} && ${RAW_MITO_DALIGNER_BIAS} -eq 1 ]]
    then
        MITO_DALIGNER_OPT="${MITO_DALIGNER_OPT} -b"
    fi
    if [[ -n ${RAW_MITO_DALIGNER_RUNID} ]]
    then
        MITO_DALIGNER_OPT="${MITO_DALIGNER_OPT} -r ${RAW_MITO_DALIGNER_RUNID}"
	else
		RAW_MITO_DALIGNER_RUNID=1
		MITO_DALIGNER_OPT="${MITO_DALIGNER_OPT} -r ${RAW_MITO_DALIGNER_RUNID}"
    fi
    if [[ -n ${RAW_MITO_DALIGNER_OLEN} ]]
    then
        MITO_DALIGNER_OPT="${MITO_DALIGNER_OPT} -l ${RAW_MITO_DALIGNER_OLEN}"
    fi    
    if [[ -n ${RAW_MITO_DALIGNER_MEM} && ${RAW_MITO_DALIGNER_MEM} -gt 0 ]]
    then
        MITO_DALIGNER_OPT="${MITO_DALIGNER_OPT} -M ${RAW_MITO_DALIGNER_MEM}"
    fi    
    if [[ -n ${RAW_MITO_DALIGNER_HITS} ]]
    then
        MITO_DALIGNER_OPT="${MITO_DALIGNER_OPT} -h ${RAW_MITO_DALIGNER_HITS}"
    fi        
    if [[ -n ${RAW_MITO_DALIGNER_T} ]]
    then
        MITO_DALIGNER_OPT="${MITO_DALIGNER_OPT} -t ${RAW_MITO_DALIGNER_T}"
    fi  
    if [[ -n ${RAW_MITO_DALIGNER_MASK} ]]
    then
        for x in ${RAW_MITO_DALIGNER_MASK}
        do 
            MITO_DALIGNER_OPT="${MITO_DALIGNER_OPT} -m ${x}"
        done
    fi
    if [[ -n ${RAW_MITO_DALIGNER_TRACESPACE} && ${RAW_MITO_DALIGNER_TRACESPACE} -gt 0 ]]
    then
        MITO_DALIGNER_OPT="${MITO_DALIGNER_OPT} -s ${RAW_MITO_DALIGNER_TRACESPACE}"
    fi
    if [[ -n ${THREADS_daligner} ]]
    then 
        MITO_DALIGNER_OPT="${MITO_DALIGNER_OPT} -j ${THREADS_daligner}"
    fi
}

myTypes=("1-mitoPrepareInput, 2-daligner, 3-LAmerge")
if [[ ${RAW_MITO_TYPE} -eq 0 ]]
then 
    ### 1-mitoPrepareInput
    if [[ ${currentStep} -eq 1 ]]
    then
        ### clean up plans 
        for x in $(ls mito_01_*_*_${RAW_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
		if [[ ! -f "${RAW_MITO_REFFASTA}" ]]
        then
        	(>&2 echo "ERROR - set RAW_MITO_REFFASTA to reference mitochondrium fasta file")
        	exit 1
   		fi        
        
        echo "${MARVEL_PATH}/bin/FA2db -v -a ${RAW_DB%.db}.db ${RAW_MITO_REFFASTA}" > mito_01_mitoPrepareInput_single_${RAW_DB%.db}.${slurmID}.plan 
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mito_01_mitoPrepareInput_single_${RAW_DB%.db}.${slurmID}.version
    ### 2-daligner
    elif [[ ${currentStep} -eq 2 ]]
    then
		### clean up plans 
        for x in $(ls mito_02_*_*_${RAW_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        ### todo get this from a config file or derive it directly from the system
        if [[ -n ${RAW_MITO_DALIGNER_NUMACTL} && ${RAW_MITO_DALIGNER_NUMACTL} -gt 0 ]]
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
            
        setDalignerOptions
        
        for x in $(seq 1 $((${rawblocks}-1)))
        do
        	echo "${NUMACTL}${MARVEL_PATH}/bin/daligner${MITO_DALIGNER_OPT} -A ${RAW_DB%.db}.${rawblocks} ${RAW_DB%.db}.${x}"        	
		done > mito_02_daligner_block_${RAW_DB%.db}.${slurmID}.plan
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mito_02_daligner_block_${RAW_DB%.db}.${slurmID}.version
	### 3-LAmerge
    elif [[ ${currentStep} -eq 3 ]]
    then
		### clean up plans 
        for x in $(ls mito_03_*_*_${RAW_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
		### call daligner options to ensure that variable RAW_MITO_DALIGNER_RUNID is set (even its not specified in the config file)
        setDalignerOptions
        
        ### create LAmerge commands 
    	echo "${MARVEL_PATH}/bin/LAmerge -n 32 ${RAW_DB%.db} ${RAW_DB%.db}.${rawblocks}.mito.las $(getSubDirName ${RAW_MITO_DALIGNER_RUNID} ${rawblocks})" > mito_03_LAmerge_single_${RAW_DB%.db}.${slurmID}.plan
    	      
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mito_03_LAmerge_single_${RAW_DB%.db}.${slurmID}.version                      	    	     
    else
        (>&2 echo "step ${currentStep} in RAW_MITO_TYPE ${RAW_MITO_TYPE} not supported")
        (>&2 echo "valid steps are: ${myTypes[${RAW_MITO_TYPE}]}")
        exit 1            
    fi    		
else
    (>&2 echo "unknown RAW_MITO_TYPE ${RAW_MITO_TYPE}")
    (>&2 echo "supported types")
    x=0; while [ $x -lt ${#myTypes[*]} ]; do (>&2 echo "${myTypes[${x}]}"); done 
    exit 1
fi

exit 0