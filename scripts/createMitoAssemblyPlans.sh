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

function setDalignerOptions()
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

function setLAfilterMitoOptions()
{
	MITO_LAFILTERMITO_OPT=" -p"		### enable purge by default
    if [[ -n ${RAW_MITO_LAFILTERMITO_VERBOSE} && ${RAW_MITO_LAFILTERMITO_VERBOSE} -gt 0 ]]
    then
        MITO_LAFILTERMITO_OPT="${MITO_LAFILTERMITO_OPT} -v"
    fi
    
    if [[ -n ${RAW_MITO_LAFILTERMITO_MINRLEN} && ${RAW_MITO_LAFILTERMITO_MINRLEN} -gt 0 ]]
    then
        MITO_LAFILTERMITO_OPT="${MITO_LAFILTERMITO_OPT} -l ${RAW_MITO_LAFILTERMITO_MINRLEN}"
	else
		RAW_MITO_LAFILTERMITO_MINRLEN=0
    fi
    
    if [[ -n ${RAW_MITO_LAFILTERMITO_MAXRLEN} && ${RAW_MITO_LAFILTERMITO_MAXRLEN} -gt 0 ]]
    then
        MITO_LAFILTERMITO_OPT="${MITO_LAFILTERMITO_OPT} -L ${RAW_MITO_LAFILTERMITO_MAXRLEN}"
	else
		(>&2 echo "WARNING - no maximum read length specified. Set max read length to (reference read length - 1000), to avoid fetching of reads that may have missed adapters.")
		RAW_MITO_LAFILTERMITO_MAXRLEN=$(($(grep -v -e ">" ${RAW_MITO_REFFASTA} | tr -d "\n" | wc -m)-1000))
		MITO_LAFILTERMITO_OPT="${MITO_LAFILTERMITO_OPT} -L ${RAW_MITO_LAFILTERMITO_MAXRLEN}"			
    fi
    
    if [[ ${RAW_MITO_LAFILTERMITO_MINRLEN} -gt ${RAW_MITO_LAFILTERMITO_MAXRLEN} ]]
    then 
    	(>&2 echo "ERROR - RAW_MITO_LAFILTERMITO_MINRLEN(${RAW_MITO_LAFILTERMITO_MINRLEN}) > RAW_MITO_LAFILTERMITO_MAXRLEN(${RAW_MITO_LAFILTERMITO_MAXRLEN})!")
    	exit 1
	fi
	
	if [[ -n ${RAW_MITO_LAFILTERMITO_UTIPS} && ${RAW_MITO_LAFILTERMITO_UTIPS} -gt 0 ]]
    then
        MITO_LAFILTERMITO_OPT="${MITO_LAFILTERMITO_OPT} -u ${RAW_MITO_LAFILTERMITO_UTIPS}"
    fi
	
	if [[ -n ${RAW_MITO_LAFILTERMITO_MAXGAPLEN} && ${RAW_MITO_LAFILTERMITO_MAXGAPLEN} -gt 0 ]]
    then
        MITO_LAFILTERMITO_OPT="${MITO_LAFILTERMITO_OPT} -g ${RAW_MITO_LAFILTERMITO_MAXGAPLEN}"
    fi

	if [[ -n ${RAW_MITO_LAFILTERMITO_MAXOVH} && ${RAW_MITO_LAFILTERMITO_MAXOVH} -gt 0 ]]
    then
        MITO_LAFILTERMITO_OPT="${MITO_LAFILTERMITO_OPT} -o ${RAW_MITO_LAFILTERMITO_MAXOVH}"
    fi
	
	if [[ -n ${RAW_MITO_LAFILTERMITO_PERCCOVLEN} && ${RAW_MITO_LAFILTERMITO_PERCCOVLEN} -gt 0 ]]
    then
        MITO_LAFILTERMITO_OPT="${MITO_LAFILTERMITO_OPT} -c ${RAW_MITO_LAFILTERMITO_PERCCOVLEN}"
    fi   
}

function setLAqOptions()
{
    MITO_LAQ_OPT=""
    adaptQTRIMCUTOFF=""    

    if [[ -n ${RAW_MITO_LAQ_MINSEG} && ${RAW_MITO_LAQ_MINSEG} -ne 0 ]]
    then
        MITO_LAQ_OPT="${MITO_LAQ_OPT} -s ${RAW_MITO_LAQ_MINSEG}"
    else 
        RAW_MITO_LAQ_MINSEG=25
        MITO_LAQ_OPT="${MITO_LAQ_OPT} -s ${RAW_MITO_LAQ_MINSEG}"
    fi

    if [[ -n ${RAW_MITO_LAQ_QTRIMCUTOFF} && ${RAW_MITO_LAQ_QTRIMCUTOFF} -ne 0 ]]
    then
        if [[ -n ${RAW_MITO_DALIGNER_TRACESPACE} && ${RAW_MITO_DALIGNER_TRACESPACE} -ne 100 ]]
        then 
            adaptQTRIMCUTOFF=$(echo "${RAW_MITO_LAQ_QTRIMCUTOFF}*${RAW_MITO_DALIGNER_TRACESPACE}/100+1" | bc)
            MITO_LAQ_OPT="${MITO_LAQ_OPT} -d ${adaptQTRIMCUTOFF}"
        else
            adaptQTRIMCUTOFF=${RAW_MITO_LAQ_QTRIMCUTOFF}
            MITO_LAQ_OPT="${MITO_LAQ_OPT} -d ${adaptQTRIMCUTOFF}"            
        fi
    else 
        if [[ -n ${RAW_MITO_DALIGNER_TRACESPACE} && ${RAW_MITO_DALIGNER_TRACESPACE} -ne 100 ]]
        then 
            RAW_MITO_LAQ_QTRIMCUTOFF=25
            adaptQTRIMCUTOFF=$(echo "${RAW_MITO_LAQ_QTRIMCUTOFF}*${RAW_MITO_DALIGNER_TRACESPACE}/100+1" | bc)
            MITO_LAQ_OPT="${MITO_LAQ_OPT} -d ${adaptQTRIMCUTOFF}"
        else
            adaptQTRIMCUTOFF=25
            RAW_MITO_LAQ_QTRIMCUTOFF=25
            MITO_LAQ_OPT="${MITO_LAQ_OPT} -d ${adaptQTRIMCUTOFF}"            
        fi
    fi
}

function setLAfilterOptions()
{
    MITO_LAFILTER_OPT=""
    	    
    if [[ -n ${RAW_MITO_LAFILTER_VERBOSE} && ${RAW_MITO_LAFILTER_VERBOSE} -ne 0 ]]
    then
        MITO_LAFILTER_OPT="${MITO_LAFILTER_OPT} -v"
    fi
    
    if [[ -n ${RAW_MITO_LAFILTER_PURGE} && ${RAW_MITO_LAFILTER_PURGE} -ne 0 ]]
    then
        MITO_LAFILTER_OPT="${MITO_LAFILTER_OPT} -p"
    fi

    if [[ -n ${RAW_MITO_LAFILTER_OLEN} && ${RAW_MITO_LAFILTER_OLEN} -ne 0 ]]
    then
        MITO_LAFILTER_OPT="${MITO_LAFILTER_OPT} -o ${RAW_MITO_LAFILTER_OLEN}"
    fi    

    if [[ -n ${RAW_MITO_LAFILTER_RLEN} && ${RAW_MITO_LAFILTER_RLEN} -ne 0 ]]
    then
        MITO_LAFILTER_OPT="${MITO_LAFILTER_OPT} -l ${RAW_MITO_LAFILTER_RLEN}"
    fi   

    if [[ -n ${RAW_MITO_LAFILTER_DIF} && ${RAW_MITO_LAFILTER_DIF} -ne 0 ]]
    then
        MITO_LAFILTER_OPT="${MITO_LAFILTER_OPT} -d ${RAW_MITO_LAFILTER_DIF}"
    fi

    if [[ -n ${RAW_MITO_LAFILTER_UBAS} ]]
    then
        MITO_LAFILTER_OPT="${MITO_LAFILTER_OPT} -u ${RAW_MITO_LAFILTER_UBAS}"
    fi

    if [[ -n ${RAW_MITO_LAFILTER_PRELOAD} && ${RAW_MITO_LAFILTER_PRELOAD} -ne 0 ]]
    then
        MITO_LAFILTER_OPT="${MITO_LAFILTER_OPT} -L"
    fi    

    if [[ -n ${RAW_MITO_LAFILTER_MINTIPCOV} && ${RAW_MITO_LAFILTER_MINTIPCOV} -gt 0 ]]
    then
        MITO_LAFILTER_OPT="${MITO_LAFILTER_OPT} -z ${RAW_MITO_LAFILTER_MINTIPCOV}"
    fi            

    if [[ -n ${RAW_MITO_LAFILTER_MULTIMAPPER} && ${RAW_MITO_LAFILTER_MULTIMAPPER} -gt 0 ]]
    then
        if [[ ${RAW_MITO_LAFILTER_MULTIMAPPER} -eq 1 ]]
        then
            MITO_LAFILTER_OPT="${MITO_LAFILTER_OPT} -w"
        else
            MITO_LAFILTER_OPT="${MITO_LAFILTER_OPT} -W"
        fi
    fi

    if [[ -n ${RAW_MITO_LAFILTER_REMPERCWORSTALN} && ${RAW_MITO_LAFILTER_REMPERCWORSTALN} -gt 0 ]]
    then
        MITO_LAFILTER_OPT="${MITO_LAFILTER_OPT} -Z ${RAW_MITO_LAFILTER_REMPERCWORSTALN}"
    fi
                    
    if [[ -n ${RAW_MITO_LAFILTER_EXCLUDEREADS} ]]
    then
        MITO_LAFILTER_OPT="${MITO_LAFILTER_OPT} -x ${RAW_MITO_LAFILTER_EXCLUDEREADS}"
    fi    

    if [[ -n ${RAW_MITO_LAFILTER_STITCH} && ${RAW_MITO_LAFILTER_STITCH} -gt 0 ]]
    then
        if [[ -n ${RAW_MITO_LAFILTER_STITCH_AGG} && ${RAW_MITO_LAFILTER_STITCH_AGG} -gt 0 ]]
        then
            MITO_LAFILTER_OPT="${MITO_LAFILTER_OPT} -S ${RAW_MITO_LAFILTER_STITCH}"
        else
            MITO_LAFILTER_OPT="${MITO_LAFILTER_OPT} -s ${RAW_MITO_LAFILTER_STITCH}"
        fi
    fi
    
	if [[ -n ${RAW_MITO_LAFILTER_TRIM} && ${RAW_MITO_LAFILTER_TRIM} -ne 0 ]] || [[ -n ${RAW_MITO_LAFILTER_UBAS} ]]
    then
        if [[ -z ${MITO_LAQ_OPT} ]]
        then 
            setLAqOptions
        fi               
        
        MITO_LAFILTER_OPT="${MITO_LAFILTER_OPT} -t trim0_d${RAW_MITO_LAQ_QTRIMCUTOFF}_s${RAW_MITO_LAQ_MINSEG} -T" 
    fi
}

function setLAfixOptions()
{
	MITO_LAFIX_OPT=""
    if [[ -n ${RAW_MITO_LAFIX_GAP} && ${RAW_MITO_LAQ_MINSEG} -ne 0 ]]
    then
        MITO_LAFIX_OPT="${MITO_LAFIX_OPT} -g ${RAW_MITO_LAFIX_GAP}"
    fi
    if [[ -n ${RAW_MITO_LAFIX_MLEN} && ${RAW_MITO_LAFIX_MLEN} -ne 0 ]]
    then
        MITO_LAFIX_OPT="${MITO_LAFIX_OPT} -x ${RAW_MITO_LAFIX_MLEN}"
    fi
    if [[ -n ${RAW_MITO_LAFIX_LOW_COVERAGE} && ${RAW_MITO_LAFIX_LOW_COVERAGE} -ne 0 ]]
    then
        MITO_LAFIX_OPT="${MITO_LAFIX_OPT} -l"
    fi
    if [[ -n ${RAW_MITO_LAFIX_MAXCHIMERLEN} && ${RAW_MITO_LAFIX_MAXCHIMERLEN} -ne 0 ]]
    then
        MITO_LAFIX_OPT="${MITO_LAFIX_OPT} -C${RAW_MITO_LAFIX_MAXCHIMERLEN}"
    fi
    if [[ -n ${RAW_MITO_LAFIX_MINCHIMERBORDERCOV} && ${RAW_MITO_LAFIX_MINCHIMERBORDERCOV} -ne 0 ]]
    then
        MITO_LAFIX_OPT="${MITO_LAFIX_OPT} -b${RAW_MITO_LAFIX_MINCHIMERBORDERCOV}"
    fi

    if [[ -z ${FIX_LAQ_OPT} ]]
    then
        setLAqOptions
    fi
    if [[ -n ${RAW_MITO_LAFIX_AGGCHIMERDETECT} && ${RAW_MITO_LAFIX_AGGCHIMERDETECT} -ne 0 ]]
    then
        MITO_LAFIX_OPT="${MITO_LAFIX_OPT} -a"
    fi
    if [[ -n ${RAW_MITO_LAFIX_DISCARDCHIMERS} && ${RAW_MITO_LAFIX_DISCARDCHIMERS} -ne 0 ]]
    then
        MITO_LAFIX_OPT="${MITO_LAFIX_OPT} -d"
    fi
    
    MITO_LAFIX_OPT="${MITO_LAFIX_OPT} -q q0_d${RAW_MITO_LAQ_QTRIMCUTOFF}_s${RAW_MITO_LAQ_MINSEG}"

    if [[ -n ${RAW_MITO_LAFIX_TRIM} && ${RAW_MITO_LAFIX_TRIM} -ne 0 ]]
    then
        MITO_LAFIX_OPT="${MITO_LAFIX_OPT} -t trim0_d${RAW_MITO_LAQ_QTRIMCUTOFF}_s${RAW_MITO_LAQ_MINSEG}"
    fi
    
    if [[ -n ${RAW_MITO_LAFIX_FIXCHIMERS} && ${RAW_MITO_LAFIX_FIXCHIMERS} -ne 0 ]]
    then
        MITO_LAFIX_OPT="${MITO_LAFIX_OPT} -X"
    fi
    
    if [[ -n ${RAW_MITO_LAFIX_CONVERTRACKS} ]]
    then
        for x in ${RAW_MITO_LAFIX_CONVERTRACKS}
        do
            MITO_LAFIX_OPT="${MITO_LAFIX_OPT} -c $x"
        done
    fi
}

function setForcealignOptions()
{
    MITO_FORCEALIGN_OPT=""
    if [[ -n ${RAW_MITO_FORCEALIGN_PARTIAL} && ${RAW_MITO_FORCEALIGN_PARTIAL} -ne 0 ]]
    then
        MITO_FORCEALIGN_OPT="${MITO_FORCEALIGN_OPT} --partial"
    fi
    if [[ -n ${RAW_MITO_FORCEALIGN_THREADS} && ${RAW_MITO_FORCEALIGN_THREADS} -gt 0 ]]
    then 
        MITO_FORCEALIGN_OPT="${MITO_FORCEALIGN_OPT} -t${RAW_MITO_FORCEALIGN_THREADS}"
    fi 
    if [[ -n ${RAW_MITO_FORCEALIGN_MAXDIST} && ${RAW_MITO_FORCEALIGN_MAXDIST} -gt 0 ]]
    then 
        MITO_FORCEALIGN_OPT="${MITO_FORCEALIGN_OPT} --maxdist${RAW_MITO_FORCEALIGN_MAXDIST}"
    fi 
    if [[ -n ${RAW_MITO_FORCEALIGN_BORDER} && ${RAW_MITO_FORCEALIGN_BORDER} -gt 0 ]]
    then 
        MITO_FORCEALIGN_OPT="${MITO_FORCEALIGN_OPT} --border${RAW_MITO_FORCEALIGN_BORDER}"
    fi 
    if [[ -n ${RAW_MITO_FORCEALIGN_CORRELATION} ]]
    then 
        MITO_FORCEALIGN_OPT="${MITO_FORCEALIGN_OPT} --correlation${RAW_MITO_FORCEALIGN_CORRELATION}"
    fi 
}

myTypes=("1-mitoPrepareInput, 2-mitodaligner, 3-mitoLAmerge, 4-mitoLAfilterMito, 5-mitoPrepareMitoHitDB, 6-mitoHitDBdaligner 7-mitoHitDBLAq 8-mitoHitDBLAfix 09_mitoPrepareMitoHitFixDB 10_mitoHitFixDBdaligner 11_mitoHitFixDBforcealign 12_mitoHitFixDBLAmerge 13_mitoHitFixDBLAq 14_mitoHitFixDBLAgap 15_mitoHitFixDBLAq 16_mitoHitFixDBLAfilter 17_mitoHitFixDBLAcorrect 18_mitoPrepareMitoHitCorDB,19_mitoHitCorDBdaligner, 20_mitoHitCorDBLAq, 21_mitoHitCorDBLAfilter")
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
                    
        setDalignerOptions
        
        #cleanup previous runs
        timeStamp=$(date '+%Y-%m-%d_%H-%M-%S')
        for x in d000_?????
        do
        	if [[ -d ${x} ]]
        	then
        		mv ${x} ${timeStamp}_${x}
        	fi	
    	done        
        
        for x in $(seq 1 $((${rawblocks}-1)))
        do
        	### todo get this from a config file or derive it directly from the system 
        	### for now this is valid the MPI batch and gpu partition
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
	        	
	        ## by default run in asymmetric mode and run_id 0 
        	echo "${NUMACTL}${MARVEL_PATH}/bin/daligner${MITO_DALIGNER_OPT} -A -r 0 ${RAW_DB%.db}.${rawblocks} ${RAW_DB%.db}.${x}"        	
		done > mito_02_mitodaligner_block_${RAW_DB%.db}.${slurmID}.plan
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mito_02_mitodaligner_block_${RAW_DB%.db}.${slurmID}.version
	### 3-LAmerge
    elif [[ ${currentStep} -eq 3 ]]
    then
		### clean up plans 
        for x in $(ls mito_03_*_*_${RAW_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
		### create LAmerge commands 
    	echo "${MARVEL_PATH}/bin/LAmerge -n 32 ${RAW_DB%.db} ${RAW_DB%.db}.${rawblocks}.mito.las $(getSubDirName 0 ${rawblocks})" > mito_03_mitoLAmerge_single_${RAW_DB%.db}.${slurmID}.plan
    	echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mito_03_mitoLAmerge_single_${RAW_DB%.db}.${slurmID}.version                      	    	     
    ### 4-mitoLAfilterMito
    elif [[ ${currentStep} -eq 4 ]]
    then
		### clean up plans 
        for x in $(ls mito_04_*_*_${RAW_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
    	setLAfilterMitoOptions
        
    ### create LAfilterMito commands 
    	echo "${MARVEL_PATH}/bin/LAfilterMito${MITO_LAFILTERMITO_OPT} ${RAW_DB%.db} ${RAW_DB%.db}.${rawblocks}.mito.las ${RAW_DB%.db}.${rawblocks}.mitoHits.las" > mito_04_mitoLAfilterMito_single_${RAW_DB%.db}.${slurmID}.plan    	      
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mito_04_mitoLAfilterMito_single_${RAW_DB%.db}.${slurmID}.version
    ### 5-mitoPrepareMitoHitDB 
	elif [[ ${currentStep} -eq 5 ]]
    then
		### clean up plans 
        for x in $(ls mito_05_*_*_${RAW_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        ### cleanup previous run if available
        timeStamp=$(date '+%Y-%m-%d_%H-%M-%S')
        if [[ -f ${RAW_DB%.db}.${rawblocks}.mitoHits.readids ]]
        then
        	mv ${RAW_DB%.db}.${rawblocks}.mitoHits.readids ${timeStamp}_${RAW_DB%.db}.${rawblocks}.mitoHits.readids
    	fi
    	if [[ -f ${RAW_DB%.db}.${rawblocks}.mitoHits.fasta ]]
        then
        	mv ${RAW_DB%.db}.${rawblocks}.mitoHits.fasta ${timeStamp}_${RAW_DB%.db}.${rawblocks}.mitoHits.fasta	
    	fi
    	if [[ -f ${PROJECT_ID}_MITO.db ]]
    	then
    		mv ${PROJECT_ID}_MITO.db ${timeStamp}_${PROJECT_ID}_MITO.db
    		for x in .${PROJECT_ID}_MITO.*
    		do
    			if [[ -f ${x} ]]
    			then
    				mv ${x} ${timeStamp}_${x}
    			fi
    		done
    	fi    	
                
        ### pull out read IDs
		echo "${MARVEL_PATH}/bin/LAshow -r ${RAW_DB%.db} ${RAW_DB%.db}.${rawblocks}.mitoHits.las | awk '{print \$2}' | sort -n -u > ${RAW_DB%.db}.${rawblocks}.mitoHits.readids" > mito_05_mitoPrepareMitoHitDB_single_${RAW_DB%.db}.${slurmID}.plan
		echo "${MARVEL_PATH}/bin/DBshow ${RAW_DB%.db} ${RAW_DB%.db}.${rawblocks}.mitoHits.readids > ${RAW_DB%.db}.${rawblocks}.mitoHits.fasta" >> mito_05_mitoPrepareMitoHitDB_single_${RAW_DB%.db}.${slurmID}.plan
		echo "${MARVEL_PATH}/bin/FA2db -v -x0 ${PROJECT_ID}_MITO_M ${RAW_DB%.db}.${rawblocks}.mitoHits.fasta" >> mito_05_mitoPrepareMitoHitDB_single_${RAW_DB%.db}.${slurmID}.plan		        
        
    	echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mito_04_mitoLAfilterMito_single_${RAW_DB%.db}.${slurmID}.version
    ### 6-mitoHitDBdaligner 7-mitoHitDBLAq 8-mitoHitDBLAfix 09_mitoPrepareMitoHitFixDB 10_mitoHitFixDBdaligner 11_mitoHitFixDBforcealign 12_mitoHitFixDBLAmerge 13_mitoHitFixDBLAq 14_mitoHitFixDBLAgap 15_mitoHitFixDBLAq 16_mitoHitFixDBLAfilter 17_mitoHitFixDBLAcorrect 18_mitoPrepareMitoHitCorDB,19_mitoHitCorDBdaligner, 20_mitoHitCorDBLAq, 21_mitoHitCorDBLAfilter
    elif [[ ${currentStep} -eq 6 ]]
    then
		### clean up plans 
        for x in $(ls mito_06_*_*_${RAW_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done   
        
        setDalignerOptions
                        
    	### todo get this from a config file or derive it directly from the system 
    	### for now this is valid the MPI batch and gpu partition
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

       	echo "${NUMACTL}${MARVEL_PATH}/bin/daligner${MITO_DALIGNER_OPT} ${PROJECT_ID}_MITO_M ${PROJECT_ID}_MITO_M" > mito_06_mitoHitDBdaligner_block_${RAW_DB%.db}.${slurmID}.plan
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mito_06_mitoHitDBdaligner_block_${RAW_DB%.db}.${slurmID}.version
    ### 7-mitoHitDBLAq 
    elif [[ ${currentStep} -eq 7 ]]
    then    
        ### clean up plans 
        for x in $(ls mito_07_*_*_${RAW_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        ### find and set LAq options 
        setLAqOptions
        
        echo "${MARVEL_PATH}/bin/LAq${MITO_LAQ_OPT} -T trim0_d${RAW_MITO_LAQ_QTRIMCUTOFF}_s${RAW_MITO_LAQ_MINSEG} -Q q0_d${RAW_MITO_LAQ_QTRIMCUTOFF}_s${RAW_MITO_LAQ_MINSEG} ${PROJECT_ID}_MITO_M ${PROJECT_ID}_MITO_M.las"  > mito_07_mitoHitDBLAq_single_${RAW_DB%.db}.${slurmID}.plan
    	echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mito_07_mitoHitDBLAq_single_${RAW_DB%.db}.${slurmID}.version
	### 8-mitoHitDBLAfix 
    elif [[ ${currentStep} -eq 8 ]]
    then    
        ### clean up plans 
        for x in $(ls mito_08_*_*_${RAW_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
    	setLAfixOptions
    	
    	echo "${MARVEL_PATH}/bin/LAfix${MITO_LAFIX_OPT} ${PROJECT_ID}_MITO_M ${PROJECT_ID}_MITO_M.las ${PROJECT_ID}_MITO_FIX_M.fasta" > mito_08_mitoHitDBLAfix_single_${RAW_DB%.db}.${slurmID}.plan    	
    	echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mito_08_mitoHitDBLAfix_single_${RAW_DB%.db}.${slurmID}.version
	### 09_mitoPrepareMitoHitFixDB 
    elif [[ ${currentStep} -eq 9 ]]
    then    
        ### clean up plans 
        for x in $(ls mito_09_*_*_${RAW_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        ## sanity check 
        if [[ ! -f ${PROJECT_ID}_MITO_FIX_M.fasta ]]
        then
        	(>&2 echo "Patched mito reads not available: ${PROJECT_ID}_MITO_FIX_M.fasta")
        	exit 1
    	fi
               
        ### cleanup previous run if available
        timeStamp=$(date '+%Y-%m-%d_%H-%M-%S')
    	if [[ -f ${PROJECT_ID}_MITO_FIX_M.db ]]
    	then
    		mv ${PROJECT_ID}_MITO_FIX_M.db ${timeStamp}_${PROJECT_ID}_MITO_FIX_M.db
    		for x in .${PROJECT_ID}_MITO_FIX_M.*
    		do
    			if [[ -f ${x} ]]
    			then
    				mv ${x} ${timeStamp}_${x}
    			fi
    		done
    	fi    	
                
        ### pull out read IDs
		echo "${MARVEL_PATH}/bin/FA2db -v -x0 -c source ${PROJECT_ID}_MITO_FIX_M ${PROJECT_ID}_MITO_FIX_M.fasta" > mito_09_mitoPrepareMitoHitFixDB_single_${RAW_DB%.db}.${slurmID}.plan
		echo "${DACCORD_PATH}/bin/fastaidrename < ${PROJECT_ID}_MITO_FIX_M.fasta | awk '{print \$1}' > ${PROJECT_ID}_MITO_FIX_D.fasta" >> mito_09_mitoPrepareMitoHitFixDB_single_${RAW_DB%.db}.${slurmID}.plan            
		echo "${DAZZLER_PATH}/bin/fasta2DB -v ${PROJECT_ID}_MITO_FIX_D ${PROJECT_ID}_MITO_FIX_D.fasta" >> mito_09_mitoPrepareMitoHitFixDB_single_${RAW_DB%.db}.${slurmID}.plan		                
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mito_09_mitoPrepareMitoHitFixDB_single_${RAW_DB%.db}.${slurmID}.version
        echo "DAZZLER $(git --git-dir=${DAZZLER_SOURCE_PATH}/DAZZ_DB/.git rev-parse --short HEAD)" >> mito_09_mitoPrepareMitoHitFixDB_single_${RAW_DB%.db}.${slurmID}.version
    	echo "fastaidrename $(git --git-dir=${DACCORD_SOURCE_PATH}/.git rev-parse --short HEAD)" >> mito_09_mitoPrepareMitoHitFixDB_single_${RAW_DB%.db}.${slurmID}.version
	### 10_mitoHitFixDBdaligner 
    elif [[ ${currentStep} -eq 10 ]]
    then    
        ### clean up plans 
        for x in $(ls mito_10_*_*_${RAW_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        setDalignerOptions               
                
        echo "${NUMACTL}${MARVEL_PATH}/bin/daligner${MITO_DALIGNER_OPT} ${PROJECT_ID}_MITO_FIX_M ${PROJECT_ID}_MITO_FIX_M" > mito_10_mitoHitFixDBdaligner_single_${RAW_DB%.db}.${slurmID}.plan
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mito_10_mitoHitFixDBdaligner_single_${RAW_DB%.db}.${slurmID}.version
	### 11_mitoHitFixDBforcealign 
    elif [[ ${currentStep} -eq 11 ]]
    then    
        ### clean up plans 
        for x in $(ls mito_11_*_*_${RAW_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        setForcealignOptions        
                
        echo "${DACCORD_PATH}/bin/${MITO_FORCEALIGN_OPT} ${PROJECT_ID}_MITO_FIX_D ${PROJECT_ID}_MITO_FIX_M ${PROJECT_ID}_MITO_FIX_M.las" > mito_11_mitoHitFixDBforcealign_single_${RAW_DB%.db}.${slurmID}.plan
		echo "forcealign $(git --git-dir=${DACCORD_SOURCE_PATH}/.git rev-parse --short HEAD)" > mito_11_mitoHitFixDBforcealign_single_${RAW_DB%.db}.${slurmID}.version	
	### 12_mitoHitFixDBLAmerge 
    elif [[ ${currentStep} -eq 12 ]]
    then    
        ### clean up plans 
        for x in $(ls mito_12_*_*_${RAW_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        
        ## sanity check
        if [[ ! -f ${PROJECT_ID}_MITO_FIX_M_f.las || ! ${PROJECT_ID}_MITO_FIX_M_r.las ]]
        then
        (>&2 echo "Missing forcealign output files: ${PROJECT_ID}_MITO_FIX_M_f.las ${PROJECT_ID}_MITO_FIX_M_r.las")
        	exit 1	
    	fi                             
                
    	echo "${MARVEL_PATH}/bin/LAmerge ${PROJECT_ID}_MITO_FIX_M ${PROJECT_ID}_MITO_FIX_M.forcealign.las ${PROJECT_ID}_MITO_FIX_M_f.las ${PROJECT_ID}_MITO_FIX_M_r.las" > mito_12_mitoHitFixDBLAmerge_single_${RAW_DB%.db}.${slurmID}.plan
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mito_12_mitoHitFixLAmerge_single_${RAW_DB%.db}.${slurmID}.version
    ### 13_mitoHitFixDBLAq
    elif [[ ${currentStep} -eq 13 ]]
    then    
        ### clean up plans 
        for x in $(ls mito_13_*_*_${RAW_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        ### find and set LAq options 
        setLAqOptions
        
        echo "${MARVEL_PATH}/bin/LAq${MITO_LAQ_OPT} -T trim0_d${RAW_MITO_LAQ_QTRIMCUTOFF}_s${RAW_MITO_LAQ_MINSEG} -Q q0_d${RAW_MITO_LAQ_QTRIMCUTOFF}_s${RAW_MITO_LAQ_MINSEG} ${PROJECT_ID}_MITO_FIX_M ${PROJECT_ID}_MITO_FIX_M.forcealign.las"  > mito_13_mitoHitFixDBLAq_single_${RAW_DB%.db}.${slurmID}.plan
    	echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mito_13_mitoHitFixDBLAq_single_${RAW_DB%.db}.${slurmID}.version        
    ### 14_mitoHitFixDBLAgap
    elif [[ ${currentStep} -eq 14 ]]
    then    
        ### clean up plans 
        for x in $(ls mito_14_*_*_${RAW_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        ### set LAq options 
        setLAqOptions
        
        echo "${MARVEL_PATH}/bin/LAgap -t trim0_d${RAW_MITO_LAQ_QTRIMCUTOFF}_s${RAW_MITO_LAQ_MINSEG} ${PROJECT_ID}_MITO_FIX_M ${PROJECT_ID}_MITO_FIX_M.forcealign.las ${PROJECT_ID}_MITO_FIX_M.gap.las" > mito_14_mitoHitFixDBLAgap_single_${RAW_DB%.db}.${slurmID}.plan
    	echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mito_14_mitoHitFixDBLAgap_single_${RAW_DB%.db}.${slurmID}.version        
    ### 15_mitoHitFixDBLAq 
    elif [[ ${currentStep} -eq 15 ]]
    then    
        ### clean up plans 
        for x in $(ls mito_15_*_*_${RAW_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        ### find and set LAq options 
        setLAqOptions
        
        echo "${MARVEL_PATH}/bin/LAq${MITO_LAQ_OPT} -T trim1_d${RAW_MITO_LAQ_QTRIMCUTOFF}_s${RAW_MITO_LAQ_MINSEG} -t trim0_d${RAW_MITO_LAQ_QTRIMCUTOFF}_s${RAW_MITO_LAQ_MINSEG} -q q0_d${RAW_MITO_LAQ_QTRIMCUTOFF}_s${RAW_MITO_LAQ_MINSEG} ${PROJECT_ID}_MITO_FIX_M ${PROJECT_ID}_MITO_FIX_M.gap.las" > mito_15_mitoHitFixDBLAq_single_${RAW_DB%.db}.${slurmID}.plan
    	echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mito_15_mitoHitFixDBLAq_single_${RAW_DB%.db}.${slurmID}.version        
	### 16_mitoHitFixDBLAfilter 
    elif [[ ${currentStep} -eq 16 ]]
    then    
        ### clean up plans 
        for x in $(ls mito_16_*_*_${RAW_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        ### find and set LAq options 
    	setLAfilterOptions
        
		echo "${MARVEL_PATH}/bin/LAfilter${MITO_LAFILTER_OPT} -T -t trim1_d${RAW_MITO_LAQ_QTRIMCUTOFF}_s${RAW_MITO_LAQ_MINSEG} ${PROJECT_ID}_MITO_FIX_M ${PROJECT_ID}_MITO_FIX_M.gap.las ${PROJECT_ID}_MITO_FIX_M.filt.las" > mito_16_mitoHitFixDBLAfilter_single_${RAW_DB%.db}.${slurmID}.plan
    	echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mito_16_mitoHitFixDBLAfilter_single_${RAW_DB%.db}.${slurmID}.version        
    ### 17_mitoHitFixDBLAcorrect 
    elif [[ ${currentStep} -eq 17 ]]
    then    
        ### clean up plans 
        for x in $(ls mito_17_*_*_${RAW_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        setLAqOptions
        
        echo "${MARVEL_PATH}/bin/LAcorrect -v -j1 -r ${PROJECT_ID}_MITO_FIX_M.filt.readIds.txt -q q0_d${RAW_MITO_LAQ_QTRIMCUTOFF}_s${RAW_MITO_LAQ_MINSEG} ${PROJECT_ID}_MITO_FIX_M ${PROJECT_ID}_MITO_FIX_M.filt.las ${PROJECT_ID}_MITO_COR_M" > mito_17_mitoHitFixDBLAcorrect_single_${RAW_DB%.db}.${slurmID}.plan
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mito_17_mitoHitFixDBLAcorrect_single_${RAW_DB%.db}.${slurmID}.version
    ### 18_mitoPrepareMitoHitCorDB
    elif [[ ${currentStep} -eq 18 ]]
    then    
        ### clean up plans 
        for x in $(ls mito_18_*_*_${RAW_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        ## sanity check 
        if [[ ! -f ${PROJECT_ID}_MITO_COR_M.00.fasta ]]
        then
        	(>&2 echo "Corrected mito reads not available: ${PROJECT_ID}_MITO.corrected.00.fasta")
        	exit 1
    	fi        
        
        ### cleanup previous run if available
        timeStamp=$(date '+%Y-%m-%d_%H-%M-%S')
    	if [[ -f ${PROJECT_ID}_MITO_COR.db ]]
    	then
    		mv ${PROJECT_ID}_MITO_COR.db ${timeStamp}_${PROJECT_ID}_MITO_COR.db
    		for x in .${PROJECT_ID}_MITO_COR.*
    		do
    			if [[ -f ${x} ]]
    			then
    				mv ${x} ${timeStamp}_${x}
    			fi
    		done
    	fi    	
                
        ### pull out read IDs
		echo "${MARVEL_PATH}/bin/FA2db -v -x0 -c source -c correctionq -c postrace ${PROJECT_ID}_MITO_COR_M ${PROJECT_ID}_MITO_COR_M.00.fasta" > mito_18_mitoPrepareMitoHitCorDB_single_${RAW_DB%.db}.${slurmID}.plan
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mito_18_mitoPrepareMitoHitCorDB_single_${RAW_DB%.db}.${slurmID}.version
    ### 19_mitoHitCorDBdaligner
    elif [[ ${currentStep} -eq 19 ]]
    then    
        ### clean up plans 
        for x in $(ls mito_19_*_*_${RAW_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        setDalignerOptions
        
    	### todo get this from a config file or derive it directly from the system 
    	### for now this is valid the MPI batch and gpu partition
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
	        	
       	echo "${NUMACTL}${MARVEL_PATH}/bin/daligner${MITO_DALIGNER_OPT} ${PROJECT_ID}_MITO_COR_M ${PROJECT_ID}_MITO_COR_M" > mito_19_mitoHitCorDBdaligner_single_${RAW_DB%.db}.${slurmID}.plan
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mito_19_mitoHitCorDBdaligner_single_${RAW_DB%.db}.${slurmID}.version
    ### 20_mitoHitCorDBLAq
    elif [[ ${currentStep} -eq 20 ]]
    then    
        ### clean up plans 
        for x in $(ls mito_20_*_*_${RAW_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        ### find and set LAq options 
        setLAqOptions
        
        echo "${MARVEL_PATH}/bin/LAq${MITO_LAQ_OPT} -T trim0_d${RAW_MITO_LAQ_QTRIMCUTOFF}_s${RAW_MITO_LAQ_MINSEG} -Q q0_d${RAW_MITO_LAQ_QTRIMCUTOFF}_s${RAW_MITO_LAQ_MINSEG} ${PROJECT_ID}_MITO_COR_M ${PROJECT_ID}_MITO_COR_M.las" > mito_20_mitoHitCorDBLAq_single_${RAW_DB%.db}.${slurmID}.plan
    	echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mito_20_mitoHitCorDBLAq_single_${RAW_DB%.db}.${slurmID}.version
    ### 21_mitoHitCorDBLAfilter
    elif [[ ${currentStep} -eq 21 ]]
    then    
        ### clean up plans 
        for x in $(ls mito_21_*_*_${RAW_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
    	setLAfilterOptions
    	
    	echo "${MARVEL_PATH}/bin/LAfilter${MITO_LAFILTER_OPT} ${PROJECT_ID}_MITO_COR_M ${PROJECT_ID}_MITO_COR_M.las ${PROJECT_ID}_MITO_COR_M.filt.las" > mito_21_mitoHitCorDBLAfilter_single_${RAW_DB%.db}.${slurmID}.plan
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mito_21_mitoHitCorDBLAfilter_single_${RAW_DB%.db}.${slurmID}.version
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