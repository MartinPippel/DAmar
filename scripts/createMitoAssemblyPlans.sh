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
    if [[ "x$1" == "x1" ]]  # use different values for polished reads !!!!!
    then 
        MITO_DALIGNER_OPT="-k 25 -e 0.92" 
    else
        if [[ -n ${RAW_MITO_DALIGNER_KMER} && ${RAW_MITO_DALIGNER_KMER} -gt 0 ]]
        then
            MITO_DALIGNER_OPT="${MITO_DALIGNER_OPT} -k ${RAW_MITO_DALIGNER_KMER}"
        fi
        if [[ -n ${RAW_MITO_DALIGNER_ERR} ]]
        then
            MITO_DALIGNER_OPT="${MITO_DALIGNER_OPT} -e ${RAW_MITO_DALIGNER_ERR}"
        fi
    fi
    
    if [[ -n ${RAW_MITO_DALIGNER_IDENTITY_OVLS} && ${RAW_MITO_DALIGNER_IDENTITY_OVLS} -gt 0 ]]
    then
        MITO_DALIGNER_OPT="${MITO_DALIGNER_OPT} -I"
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
    
    if [[ "x$1" == "x1" ]]
    then
        ### for corercted reads allow zero-qvalued aliugn segments
        MITO_LAQ_OPT="${MITO_LAQ_OPT} -c"
		if [[ -n ${COR_MITO_LAQ_MINSEG} && ${COR_MITO_LAQ_MINSEG} -ne 0 ]]
	    then
	        MITO_LAQ_OPT="${COR_LAQ_OPT} -s ${COR_MITO_LAQ_MINSEG}"
	    else 
	        COR_MITO_LAQ_MINSEG=2
	        MITO_LAQ_OPT="${MITO_LAQ_OPT} -s ${COR_MITO_LAQ_MINSEG}"
	    fi
	
	    if [[ -n ${COR_MITO_LAQ_QTRIMCUTOFF} && ${COR_MITO_LAQ_QTRIMCUTOFF} -ne 0 ]]
	    then
	        if [[ -n ${RAW_MITO_DALIGNER_TRACESPACE} && ${RAW_MITO_DALIGNER_TRACESPACE} -ne 100 ]]
	        then 
	            adaptQTRIMCUTOFF=$(echo "${COR_MITO_LAQ_QTRIMCUTOFF}*${RAW_MITO_DALIGNER_TRACESPACE}/100+1" | bc)
	            MITO_LAQ_OPT="${MITO_LAQ_OPT} -d ${adaptQTRIMCUTOFF}"
	        else
	            adaptQTRIMCUTOFF=${COR_MITO_LAQ_QTRIMCUTOFF}
	            MITO_LAQ_OPT="${MITO_LAQ_OPT} -d ${adaptQTRIMCUTOFF}"            
	        fi
	    else 
	        if [[ -n ${COR_MITO_DALIGNER_TRACESPACE} && ${COR_MITO_DALIGNER_TRACESPACE} -ne 100 ]]
	        then 
	            COR_MITO_LAQ_QTRIMCUTOFF=25
	            adaptQTRIMCUTOFF=$(echo "${COR_MITO_LAQ_QTRIMCUTOFF}*${RAW_MITO_DALIGNER_TRACESPACE}/100+1" | bc)
	            MITO_LAQ_OPT="${MITO_LAQ_OPT} -d ${adaptQTRIMCUTOFF}"
	        else
	            adaptQTRIMCUTOFF=25
	            COR_MITO_LAQ_QTRIMCUTOFF=25
	            MITO_LAQ_OPT="${MITO_LAQ_OPT} -d ${adaptQTRIMCUTOFF}"            
	        fi
	    fi    	
	else
	    if [[ -n ${RAW_MITO_LAQ_MINSEG} && ${RAW_MITO_LAQ_MINSEG} -ne 0 ]]
	    then
	        MITO_LAQ_OPT="${MITO_LAQ_OPT} -s ${RAW_MITO_LAQ_MINSEG}"
	    else 
	        RAW_MITO_LAQ_MINSEG=2
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
        if [[ "x$1" == "x1" ]]
        then 
            setLAqOptions 1
    	else
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

function setDaccordOptions()
{
	MITO_DACCORD_OPT=""
	
	if [[ -z ${RAW_MITO_DACCORD_THREADS} ]]
	then 
		RAW_MITO_DACCORD_THREADS=8
	fi
	MITO_DACCORD_OPT="${MITO_DACCORD_OPT} -t${RAW_MITO_DACCORD_THREADS}"
	
	if [[ -n ${RAW_MITO_DACCORD_WINDOW} && ${RAW_MITO_DACCORD_WINDOW} -gt 0 ]]
	then 
		MITO_DACCORD_OPT="${MITO_DACCORD_OPT} -w${RAW_MITO_DACCORD_WINDOW}"
	fi

	if [[ -n ${RAW_MITO_DACCORD_ADVANCESIZE} && ${RAW_MITO_DACCORD_ADVANCESIZE} -gt 0 ]]
	then 
		MITO_DACCORD_OPT="${MITO_DACCORD_OPT} -a${RAW_MITO_DACCORD_ADVANCESIZE}"
	fi
	
	if [[ -n ${RAW_MITO_DACCORD_MAXDEPTH} && ${RAW_MITO_DACCORD_MAXDEPTH} -gt 0 ]]
	then 
		MITO_DACCORD_OPT="${MITO_DACCORD_OPT} -d${RAW_MITO_DACCORD_MAXDEPTH}"
	fi
	
	if [[ -n ${RAW_MITO_DACCORD_FULLSEQ} && ${RAW_MITO_DACCORD_FULLSEQ} -gt 0 ]]
	then 
		MITO_DACCORD_OPT="${MITO_DACCORD_OPT} -f1"
	fi
	
	if [[ -n ${RAW_MITO_DACCORD_VEBOSE} && ${RAW_MITO_DACCORD_VEBOSE} -gt 0 ]]
	then 
		MITO_DACCORD_OPT="${MITO_DACCORD_OPT} -V${RAW_MITO_DACCORD_VEBOSE}"
	fi
		
	if [[ -n ${RAW_MITO_DACCORD_MINWINDOWCOV} && ${RAW_MITO_DACCORD_MINWINDOWCOV} -gt 0 ]]
	then 
		MITO_DACCORD_OPT="${MITO_DACCORD_OPT} -m${RAW_MITO_DACCORD_MINWINDOWCOV}"
	fi
	
	if [[ -n ${RAW_MITO_DACCORD_MINWINDOWERR} && ${RAW_MITO_DACCORD_MINWINDOWERR} -gt 0 ]]
	then 
		MITO_DACCORD_OPT="${MITO_DACCORD_OPT} -e${RAW_MITO_DACCORD_MINWINDOWERR}"
	fi
	
	if [[ -n ${RAW_MITO_DACCORD_MINOUTLEN} && ${RAW_MITO_DACCORD_MINOUTLEN} -gt 0 ]]
	then 
		MITO_DACCORD_OPT="${MITO_DACCORD_OPT} -l${RAW_MITO_DACCORD_MINOUTLEN}"
	fi
	
	if [[ -n ${RAW_MITO_DACCORD_MINKFREQ} && ${RAW_MITO_DACCORD_MINKFREQ} -gt 0 ]]
	then 
		MITO_DACCORD_OPT="${MITO_DACCORD_OPT} --minfilterfreq${RAW_MITO_DACCORD_MINKFREQ}"
	fi
	
	if [[ -n ${RAW_MITO_DACCORD_MAXKFREQ} && ${RAW_MITO_DACCORD_MAXKFREQ} -gt 0 ]]
	then 
		MITO_DACCORD_OPT="${MITO_DACCORD_OPT} --maxfilterfreq${RAW_MITO_DACCORD_MAXKFREQ}"
	fi
	
	if [[ -n ${RAW_MITO_DACCORD_MAXOVLS} && ${RAW_MITO_DACCORD_MAXOVLS} -gt 0 ]]
	then 
		MITO_DACCORD_OPT="${MITO_DACCORD_OPT} -D${RAW_MITO_DACCORD_MAXOVLS}"
	fi
	
	if [[ -n ${RAW_MITO_DACCORD_VARD} && ${RAW_MITO_DACCORD_VARD} -gt 0 ]]
	then 
		MITO_DACCORD_OPT="${MITO_DACCORD_OPT} --vard${RAW_MITO_DACCORD_VARD}"
	fi
	
	if [[ -n ${RAW_MITO_DACCORD_KMER} && ${RAW_MITO_DACCORD_KMER} -gt 0 ]]
	then 
		MITO_DACCORD_OPT="${MITO_DACCORD_OPT} -k${RAW_MITO_DACCORD_KMER}"
	fi
}

# type-0 steps [1-21]: 1-mitoPrepareInput, 2-mitodaligner, 3-mitoLAmerge, 4-mitoLAfilterMito, 5-mitoPrepareMitoHitDB, 6-mitoHitDBdaligner 7-mitoHitDBLAq 8-mitoHitDBLAfix 09_mitoPrepareMitoHitFixDB, 10_mitoHitFixDBdaligner, 11_mitoHitFixDBLAq, 
#                      12_mitoHitFixDBLAgap 13_mitoHitFixDBLAq 14_mitoHitFixDBLAfilter 15_mitoHitFixDBLAcorrect 16_mitoPrepareMitoHitCorDB 17_mitoHitCorDBdaligner, 18_mitoHitCorDBLAq, 19_mitoHitCorDBLAfilter, 20_mitoHitCorDBTour 21_mitoHitCorDBArrowPolishing 
myTypes=("1-mitoPrepareInput 2-mitodaligner 3-mitoLAmerge 4-mitoLAfilterMito 5-mitoPrepareMitoHitDB 6-mitoHitDBdaligner 7-mitoHitDBLAq 8-mitoHitDBLAfix 09_mitoPrepareMitoHitFixDB 10_mitoHitFixDBdaligner 11_mitoHitFixDBLAq 12_mitoHitFixDBLAgap 13_mitoHitFixDBLAq 14_mitoHitFixDBLAfilter 15_mitoHitFixDBLAcorrect 16_mitoPrepareMitoHitCorDB 17_mitoHitCorDBdaligner 18_mitoHitCorDBLAq 19_mitoHitCorDBLAfilter 20_mitoHitCorDBTour 21_mitoHitCorDBArrowPolishing")
if [[ ${RAW_MITO_TYPE} -eq 0 ]]
then
	if [[ ${currentStep} -lt 10 ]]
	then 
		sID=0${currentStep}
	else
		sID=${currentStep}
	fi
	
    ### 1-mitoPrepareInput
    if [[ ${currentStep} -eq 1 ]]
    then
        ### clean up plans
        for x in $(ls mito_${sID}_*_*_${RAW_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
		if [[ ! -f "${RAW_MITO_REFFASTA}" ]]
        then
        	(>&2 echo "ERROR - set RAW_MITO_REFFASTA to reference mitochondrium fasta file")
        	exit 1
   		fi        
        
        echo "${MARVEL_PATH}/bin/FA2db -v -a ${RAW_DB%.db}.db ${RAW_MITO_REFFASTA}" > mito_${sID}_mitoPrepareInput_single_${RAW_DB%.db}.${slurmID}.plan 
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mito_${sID}_mitoPrepareInput_single_${RAW_DB%.db}.${slurmID}.version
    ### 2-daligner
    elif [[ ${currentStep} -eq 2 ]]
    then
		### clean up plans 
        for x in $(ls mito_${sID}_*_*_${RAW_DB}.${slurmID}.* 2> /dev/null)
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
	        ## by default run in asymmetric mode and run_id 0 
        	echo "${MARVEL_PATH}/bin/daligner${MITO_DALIGNER_OPT} -A -r 0 ${RAW_DB%.db}.${rawblocks} ${RAW_DB%.db}.${x}"        	
		done > mito_${sID}_mitodaligner_block_${RAW_DB%.db}.${slurmID}.plan
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mito_${sID}_mitodaligner_block_${RAW_DB%.db}.${slurmID}.version
	### 3-LAmerge
    elif [[ ${currentStep} -eq 3 ]]
    then
		### clean up plans 
        for x in $(ls mito_${sID}_*_*_${RAW_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
		### create LAmerge commands 
    	echo "${MARVEL_PATH}/bin/LAmerge -n 32 ${RAW_DB%.db} ${RAW_DB%.db}.${rawblocks}.mito.las $(getSubDirName 0 ${rawblocks})" > mito_${sID}_mitoLAmerge_single_${RAW_DB%.db}.${slurmID}.plan
    	echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mito_${sID}_mitoLAmerge_single_${RAW_DB%.db}.${slurmID}.version                      	    	     
    ### 4-mitoLAfilterMito
    elif [[ ${currentStep} -eq 4 ]]
    then
		### clean up plans 
        for x in $(ls mito_${sID}_*_*_${RAW_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
    	setLAfilterMitoOptions
        
    ### create LAfilterMito commands 
    	echo "${MARVEL_PATH}/bin/LAfilterMito${MITO_LAFILTERMITO_OPT} ${RAW_DB%.db} ${RAW_DB%.db}.${rawblocks}.mito.las ${RAW_DB%.db}.${rawblocks}.mitoHits.las" > mito_${sID}_mitoLAfilterMito_single_${RAW_DB%.db}.${slurmID}.plan    	      
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mito_${sID}_mitoLAfilterMito_single_${RAW_DB%.db}.${slurmID}.version
    ### 5-mitoPrepareMitoHitDB 
	elif [[ ${currentStep} -eq 5 ]]
    then
		### clean up plans 
        for x in $(ls mito_${sID}_*_*_${RAW_DB}.${slurmID}.* 2> /dev/null)
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
		echo "${MARVEL_PATH}/bin/LAshow -r ${RAW_DB%.db} ${RAW_DB%.db}.${rawblocks}.mitoHits.las | awk '{print \$2}' | sort -n -u > ${RAW_DB%.db}.${rawblocks}.mitoHits.readids" > mito_${sID}_mitoPrepareMitoHitDB_single_${RAW_DB%.db}.${slurmID}.plan
		echo "${MARVEL_PATH}/bin/DBshow ${RAW_DB%.db} ${RAW_DB%.db}.${rawblocks}.mitoHits.readids > ${RAW_DB%.db}.${rawblocks}.mitoHits.fasta" >> mito_${sID}_mitoPrepareMitoHitDB_single_${RAW_DB%.db}.${slurmID}.plan
		echo "${MARVEL_PATH}/bin/FA2db -v -x0 ${PROJECT_ID}_MITO_M ${RAW_DB%.db}.${rawblocks}.mitoHits.fasta" >> mito_${sID}_mitoPrepareMitoHitDB_single_${RAW_DB%.db}.${slurmID}.plan
        echo "${MARVEL_PATH}/bin/DBsplit ${PROJECT_ID}_MITO_M" >> mito_${sID}_mitoPrepareMitoHitDB_single_${RAW_DB%.db}.${slurmID}.plan
        echo "${DACCORD_PATH}/bin/fastaidrename < ${RAW_DB%.db}.${rawblocks}.mitoHits.fasta | awk '{print \$1}' > ${RAW_DB%.db}.${rawblocks}.mitoHitsDazzler.fasta" >> mito_${sID}_mitoPrepareMitoHitDB_single_${RAW_DB%.db}.${slurmID}.plan
		echo "${DAZZLER_PATH}/bin/fasta2DB -v ${PROJECT_ID}_MITO_D ${RAW_DB%.db}.${rawblocks}.mitoHitsDazzler.fasta" >> mito_${sID}_mitoPrepareMitoHitDB_single_${RAW_DB%.db}.${slurmID}.plan
        echo "${DAZZLER_PATH}/bin/DBsplit ${PROJECT_ID}_MITO_D" >> mito_${sID}_mitoPrepareMitoHitDB_single_${RAW_DB%.db}.${slurmID}.plan

    	echo "MARVEL LAshow $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mito_${sID}_mitoPrepareMitoHitDB_single_${RAW_DB%.db}.${slurmID}.version
    elif [[ ${currentStep} -eq 6 ]]
    then
		### clean up plans 
        for x in $(ls mito_${sID}_*_*_${RAW_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done   
        
        setDalignerOptions                    

       	echo "${MARVEL_PATH}/bin/daligner${MITO_DALIGNER_OPT} ${PROJECT_ID}_MITO_M ${PROJECT_ID}_MITO_M" > mito_${sID}_mitoHitDBdaligner_block_${RAW_DB%.db}.${slurmID}.plan
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mito_${sID}_mitoHitDBdaligner_block_${RAW_DB%.db}.${slurmID}.version
    ### 07-lassort
	elif [[ ${currentStep} -eq 7 ]]
    then
        ### clean up plans 
        for x in $(ls mito_${sID}_*_*_${RAW_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done   
        		
        echo "${LASTOOLS_PATH}/bin/lassort -sfull ${PROJECT_ID}_MITO_M.sort.las ${PROJECT_ID}_MITO_M.las" > mito_${sID}_mitoLassort_single_${RAW_DB%.db}.${slurmID}.plan
        echo "LASTOOLS lassort $(git --git-dir=${DACCORD_SOURCE_PATH}/.git rev-parse --short HEAD)" > mito_${sID}_mitoLassort_single_${RAW_DB%.db}.${slurmID}.version
    ### 08-computeIntrinsicQV
	elif [[ ${currentStep} -eq 8 ]]
    then
        ### clean up plans 
        for x in $(ls mito_${sID}_*_*_${RAW_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done

        cov=$(tail -n 1 ${PROJECT_ID}_MITO_M.db | tr -d " ")
        echo "${DACCORD_PATH}/bin/computeintrinsicqv2 -d${cov} ${PROJECT_ID}_MITO_D.db ${PROJECT_ID}_MITO_M.sort.las" > mito_${sID}_mitocomputeintrinsicqv2_single_${RAW_DB%.db}.${slurmID}.plan
		echo "DACCORD computeintrinsicqv2 $(git --git-dir=${DACCORD_SOURCE_PATH}/.git rev-parse --short HEAD)" > mito_${sID}_mitocomputeintrinsicqv2_single_${RAW_DB%.db}.${slurmID}.version

    ### 09-daccord
    elif [[ ${currentStep} -eq 9 ]]
    then    
        ### clean up plans 
        for x in $(ls mito_${sID}_*_*_${RAW_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        ### find and set LAq options 
        setDaccordOptions 
    
        echo "${DACCORD_PATH}/bin/daccord ${MITO_DACCORD_OPT} --eprofonly -E${PROJECT_ID}_MITO_M.sort.eprof ${PROJECT_ID}_MITO_M.sort.las ${PROJECT_ID}_MITO_D.db && ${DACCORD_PATH}/bin/daccord ${MITO_DACCORD_OPT} -E${PROJECT_ID}_MITO_M.sort.eprof ${PROJECT_ID}_MITO_M.sort.las ${PROJECT_ID}_MITO_D.db > ${PROJECT_ID}_MITO_M.sort.dac.fasta" > mito_${sID}_mitodaccord_single_${RAW_DB%.db}.${slurmID}.plan
    	echo "DACCORD daccord $(git --git-dir=${DACCORD_SOURCE_PATH}/.git rev-parse --short HEAD)" > mito_${sID}_mitodaccord_single_${RAW_DB%.db}.${slurmID}.version
    ### 10_mitoPrepareMitoHitCorDB
    elif [[ ${currentStep} -eq 10 ]]
    then    
        ### clean up plans 
        for x in $(ls mito_${sID}_*_*_${RAW_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        ## sanity check 
        if [[ ! -f ${PROJECT_ID}_MITO_M.sort.dac.fasta ]]
        then
        	(>&2 echo "Corrected mito reads not available: ${PROJECT_ID}_MITO_M.sort.dac.fasta")
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
                
        echo "${MARVEL_PATH}/scripts/splitReads.py ${PROJECT_ID}_MITO_M.sort.dac.fasta ${PROJECT_ID}_MITO_M.sort.dac.split.fasta 5000 0" > mito_${sID}_mitoPrepareMitoHitCorDB_single_${RAW_DB%.db}.${slurmID}.plan
		echo "${MARVEL_PATH}/bin/FA2db -v -x0 ${PROJECT_ID}_MITO_COR_M ${PROJECT_ID}_MITO_M.sort.dac.split.fasta" >> mito_${sID}_mitoPrepareMitoHitCorDB_single_${RAW_DB%.db}.${slurmID}.plan
		echo "${DAZZLER_PATH}/bin/fasta2DB -v ${PROJECT_ID}_MITO_COR_D ${PROJECT_ID}_MITO_M.sort.dac.split.fasta" >> mito_${sID}_mitoPrepareMitoHitCorDB_single_${RAW_DB%.db}.${slurmID}.plan		                
		echo "${MARVEL_PATH}/bin/DBsplit ${PROJECT_ID}_MITO_COR_M" >> mito_${sID}_mitoPrepareMitoHitCorDB_single_${RAW_DB%.db}.${slurmID}.plan
		echo "${DAZZLER_PATH}/bin/DBsplit ${PROJECT_ID}_MITO_COR_D" >> mito_${sID}_mitoPrepareMitoHitCorDB_single_${RAW_DB%.db}.${slurmID}.plan		                
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mito_${sID}_mitoPrepareMitoHitCorDB_single_${RAW_DB%.db}.${slurmID}.version
        echo "DAZZLER $(git --git-dir=${DAZZLER_SOURCE_PATH}/DAZZ_DB/.git rev-parse --short HEAD)" >> mito_${sID}_mitoPrepareMitoHitCorDB_single_${RAW_DB%.db}.${slurmID}.version
    ### 11_mitoHitCorDBdaligner
    elif [[ ${currentStep} -eq 11 ]]
    then    
        ### clean up plans 
        for x in $(ls mito_${sID}_*_*_${RAW_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        setDalignerOptions 1
        
    	echo "${MARVEL_PATH}/bin/daligner${MITO_DALIGNER_OPT} ${PROJECT_ID}_MITO_COR_M ${PROJECT_ID}_MITO_COR_M" > mito_${sID}_mitoHitCorDBdaligner_single_${RAW_DB%.db}.${slurmID}.plan
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mito_${sID}_mitoHitCorDBdaligner_single_${RAW_DB%.db}.${slurmID}.version
    ### 12_mitoHitCorDBLAq
    elif [[ ${currentStep} -eq 12 ]]
    then    
        ### clean up plans 
        for x in $(ls mito_${sID}_*_*_${RAW_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        ### find and set LAq options 
        setLAqOptions 1
        
        echo "${MARVEL_PATH}/bin/LAq${MITO_LAQ_OPT} -T trim0_d${COR_MITO_LAQ_QTRIMCUTOFF}_s${COR_MITO_LAQ_MINSEG} -Q q0_d${COR_MITO_LAQ_QTRIMCUTOFF}_s${COR_MITO_LAQ_MINSEG} ${PROJECT_ID}_MITO_COR_M ${PROJECT_ID}_MITO_COR_M.las" > mito_${sID}_mitoHitCorDBLAq_single_${RAW_DB%.db}.${slurmID}.plan
    	echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mito_${sID}_mitoHitCorDBLAq_single_${RAW_DB%.db}.${slurmID}.version
    ### 13_mitoHitCorDBLAfilter 
    elif [[ ${currentStep} -eq 13 ]]
    then    
        ### clean up plans 
        for x in $(ls mito_${sID}_*_*_${RAW_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
    	setLAfilterOptions 1
    	
    	echo "${MARVEL_PATH}/bin/LAfilter${MITO_LAFILTER_OPT} ${PROJECT_ID}_MITO_COR_M ${PROJECT_ID}_MITO_COR_M.las ${PROJECT_ID}_MITO_COR_M.filt.las" > mito_${sID}_mitoHitCorDBLAfilter_single_${RAW_DB%.db}.${slurmID}.plan
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mito_${sID}_mitoHitCorDBLAfilter_single_${RAW_DB%.db}.${slurmID}.version
    ### 14_mitoHitCorDBTour     
    elif [[ ${currentStep} -eq 14 ]]
    then
    	### clean up plans 
        for x in $(ls mito_${sID}_*_*_${RAW_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done    
        
        setLAfilterOptions 1
        
        ### find and set LAq options 
        setLAqOptions 1
        
        echo "${MARVEL_PATH}/bin/OGbuild -t trim0_d${COR_MITO_LAQ_QTRIMCUTOFF}_s${COR_MITO_LAQ_MINSEG} ${PROJECT_ID}_MITO_COR_M ${PROJECT_ID}_MITO_COR_M.filt.las ${PROJECT_ID}_MITO_COR_M.graphml" > mito_${sID}_mitoHitCorDBTour_single_${RAW_DB%.db}.${slurmID}.plan
        echo "${MARVEL_PATH}/scripts/OGtour.py -c -d -l3 -b4 ${PROJECT_ID}_MITO_COR_M ${PROJECT_ID}_MITO_COR_M.graphml" >> mito_${sID}_mitoHitCorDBTour_single_${RAW_DB%.db}.${slurmID}.plan
        echo "${MARVEL_PATH}/scripts/tour2fasta.py -t trim0_d${COR_MITO_LAQ_QTRIMCUTOFF}_s${COR_MITO_LAQ_MINSEG} -p ${PROJECT_ID}_MITO_COR_M ${PROJECT_ID}_MITO_COR_M ${PROJECT_ID}_MITO_COR_M.graphml ${PROJECT_ID}_MITO_COR_M.tour.paths" >> mito_${sID}_mitoHitCorDBTour_single_${RAW_DB%.db}.${slurmID}.plan
        echo "${MARVEL_PATH}/bin/OGlayout -R -d 300 ${PROJECT_ID}_MITO_COR_M.tour.graphml ${PROJECT_ID}_MITO_COR_M.tour.layout.graphml" >> mito_${sID}_mitoHitCorDBTour_single_${RAW_DB%.db}.${slurmID}.plan
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mito_${sID}_mitoHitCorDBTour_single_${RAW_DB%.db}.${slurmID}.version
    ### 15_mitoHitCorDBArrowPolishing    
    elif [[ ${currentStep} -eq 15 ]]
    then
    	### clean up plans 
        for x in $(ls mito_${sID}_*_*_${RAW_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done    
        
        ## TODO --> pull out raw read ids, create pacbio data set, run arrow polishing twice !!!!!!
        
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

