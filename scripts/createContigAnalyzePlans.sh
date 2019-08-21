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

if [[ ! -n "${FIX_CORR_TYPE}" ]]
then 
    (>&2 echo "cannot create touring scripts if variable FIX_CORR_TYPE is not set.")
    exit 1
fi

if [[ ! -n ${RAW_DB} ]]
then 
    (>&2 echo "raw database unknown - You have to set the variable RAW_DB")
    exit 1
fi

if [[ ! -f ${RAW_DB%.db}.db ]]
then 
    (>&2 echo "raw database ${RAW_DB%.db}.db missing")
    exit 1 
fi

if [[ ! -n ${FIX_DB} ]]
then 
    (>&2 echo "patched database unknown - You have to set the variable FIX_DB")
    exit 1
fi

if [[ ! -f ${FIX_DB%.db}.db ]]
then 
    (>&2 echo "patched database ${FIX_DB%.db}.db missing")
    exit 1
fi

if [[ ! -n ${COR_DB} ]]
then 
    (>&2 echo "corrected database unknown - You have to set the variable COR_DB")
    exit 1
fi

if [[ -z ${CONT_DB} ]]
then 
    (>&2 echo "contig database unknown - You have to set the variable CONT_DB")
    exit 1
fi

if [[ -z ${CONT_DAZZ_DB} ]]
then 
    (>&2 echo "contig database unknown - You have to set the variable CONT_DAZZ_DB")
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

function setLAfilterChainsOptions()
{
    CONTIG_LAFILTERCHAINS_OPT=""

    if [[ -n ${COR_CONTIG_LAFILTERCHAINS_NREP} && ${COR_CONTIG_LAFILTERCHAINS_NREP} -ne 0 ]]
    then
        CONTIG_LAFILTERCHAINS_OPT="${CONTIG_LAFILTERCHAINS_OPT} -n ${COR_CONTIG_LAFILTERCHAINS_NREP}"
	else
		 CONTIG_LAFILTERCHAINS_OPT="${CONTIG_LAFILTERCHAINS_OPT} -n 500" 
    fi

    if [[ -n ${COR_CONTIG_LAFILTERCHAINS_PURGE} && ${COR_CONTIG_LAFILTERCHAINS_PURGE} -ne 0 ]]
    then
        CONTIG_LAFILTERCHAINS_OPT="${CONTIG_LAFILTERCHAINS_OPT} -p"
    fi

    if [[ -n ${COR_CONTIG_LAFILTERCHAINS_KEEP} ]]
    then
        CONTIG_LAFILTERCHAINS_OPT="${CONTIG_LAFILTERCHAINS_OPT} -k ${COR_CONTIG_LAFILTERCHAINS_KEEP}"
    fi
    
    if [[ -n ${COR_CONTIG_LAFILTERCHAINS_FUZZ} && ${COR_CONTIG_LAFILTERCHAINS_FUZZ} -gt 0 ]]
    then
        CONTIG_LAFILTERCHAINS_OPT="${CONTIG_LAFILTERCHAINS_OPT} -f ${COR_CONTIG_LAFILTERCHAINS_FUZZ}"
    fi
    
    if [[ -n ${COR_CONTIG_LAFILTERCHAINS_PERCENTCONTAINED} && ${COR_CONTIG_LAFILTERCHAINS_PERCENTCONTAINED} -gt 0 && ${COR_CONTIG_LAFILTERCHAINS_PERCENTCONTAINED} -le 100 ]]
    then
        CONTIG_LAFILTERCHAINS_OPT="${CONTIG_LAFILTERCHAINS_OPT} -c ${COR_CONTIG_LAFILTERCHAINS_PERCENTCONTAINED}"
    fi
    
    if [[ -n ${COR_CONTIG_LAFILTERCTCHAIN_REPEATTRACK} ]]
    then
    	CONTIG_LAFILTERCHAINS_OPT="${CONTIG_LAFILTERCHAINS_OPT} -r ${COR_CONTIG_LAFILTERCTCHAIN_REPEATTRACK}"	
	fi
}

function setcreateCorrectedContigDBOptions()
{
	if [[ -x "${ANALYZE_DIR}" ]] 
	then
		ANALYZE_DIR=analyzeContigs	
	fi 
	
}

function setDBdustOptions()
{
	if [[ -x "${ANALYZE_DIR}" ]] 
	then
		setcreateCorrectedContigDBOptions
	fi

    CONTIG_DBDUST_OPT=""
    if [[ -n ${COR_CONTIG_DBDUST_BIAS} && ${COR_CONTIG_DBDUST_BIAS} -ge 1 ]]
    then
        CONTIG_DBDUST_OPT="${CONTIG_DBDUST_OPT} -b"
    fi
}

function setCatrackOptions()
{
	if [[ -x "${ANALYZE_DIR}" ]] 
	then
		setcreateCorrectedContigDBOptions
	fi
	
    CONTIG_CATRACK_OPT=""
    if [[ -n ${COR_CONTIG_CATRACK_VERBOSE} && ${COR_CONTIG_CATRACK_VERBOSE} -ge 1 ]]
    then
        CONTIG_CATRACK_OPT="${CONTIG_CATRACK_OPT} -v"
    fi
    if [[ -n ${COR_CONTIG_CATRACK_DELETE} && ${COR_CONTIG_CATRACK_DELETE} -ge 1 ]]
    then
        CONTIG_CATRACK_OPT="${CONTIG_CATRACK_OPT} -d"
    fi
    if [[ -n ${COR_CONTIG_CATRACK_OVERWRITE} && ${COR_CONTIG_CATRACK_OVERWRITE} -ge 1 ]]
    then
        CONTIG_CATRACK_OPT="${CONTIG_CATRACK_OPT} -f"
    fi
}

function setDatanderOptions()
{
	if [[ -x "${ANALYZE_DIR}" ]] 
	then
		setcreateCorrectedContigDBOptions
	fi
	
    ### find and set datander options 
    CONTIG_DATANDER_OPT=""
    if [[ -n ${COR_CONTIG_DATANDER_THREADS} ]]
    then
        CONTIG_DATANDER_OPT="${CONTIG_DATANDER_OPT} -j ${COR_CONTIG_DATANDER_THREADS}"
    fi
    if [[ -n ${COR_CONTIG_DATANDER_MINLEN} ]]
    then
        CONTIG_DATANDER_OPT="${CONTIG_DATANDER_OPT} -l ${COR_CONTIG_DATANDER_MINLEN}"
    fi
    if [[ -n ${COR_CONTIG_DATANDER_FOLDER} ]]
    then
        CONTIG_DATANDER_OPT="${CONTIG_DATANDER_OPT} -o ${COR_CONTIG_DATANDER_FOLDER}"
    else
        COR_CONTIG_DATANDER_FOLDER="tan"
        CONTIG_DATANDER_OPT="${CONTIG_DATANDER_OPT} -o ${COR_CONTIG_DATANDER_FOLDER}"
    fi
}
                  
function setTANmaskOptions()
{
    CONTIG_TANMASK_OPT=""
    if [[ -n ${COR_CONTIG_TANMASK_VERBOSE} && ${COR_CONTIG_TANMASK_VERBOSE} -ge 1 ]]
    then
        CONTIG_TANMASK_OPT="${CONTIG_TANMASK_OPT} -v"
    fi
    if [[ -n ${COR_CONTIG_TANMASK_MINLEN} && ${COR_CONTIG_TANMASK_MINLEN} -ge 1 ]]
    then
        CONTIG_TANMASK_OPT="${CONTIG_TANMASK_OPT} -l ${COR_CONTIG_TANMASK_MINLEN}"
    fi
    if [[ -n ${COR_CONTIG_TANMASK_TRACK} ]]
    then
        CONTIG_TANMASK_OPT="${CONTIG_TANMASK_OPT} -m ${COR_CONTIG_TANMASK_TRACK}"
    fi
}

function setDaligerOptions()
{
    CONTIG_DALIGNER_OPT=""
    if [[ -n ${COR_CONTIG_DALIGNER_IDENTITY_OVLS} && ${COR_CONTIG_DALIGNER_IDENTITY_OVLS} -gt 0 ]]
    then
        CONTIG_DALIGNER_OPT="${CONTIG_DALIGNER_OPT} -I"
    fi
    if [[ -n ${COR_CONTIG_DALIGNER_KMER} && ${COR_CONTIG_DALIGNER_KMER} -gt 0 ]]
    then
        CONTIG_DALIGNER_OPT="${CONTIG_DALIGNER_OPT} -k ${COR_CONTIG_DALIGNER_KMER}"
    fi
    if [[ -n ${COR_CONTIG_DALIGNER_ERR} ]]
    then
        CONTIG_DALIGNER_OPT="${CONTIG_DALIGNER_OPT} -e ${COR_CONTIG_DALIGNER_ERR}"
    fi
    if [[ -n ${COR_CONTIG_DALIGNER_BIAS} && ${COR_CONTIG_DALIGNER_BIAS} -eq 1 ]]
    then
        CONTIG_DALIGNER_OPT="${CONTIG_DALIGNER_OPT} -b"
    fi
    if [[ -n ${COR_CONTIG_DALIGNER_RUNID} ]]
    then
        CONTIG_DALIGNER_OPT="${CONTIG_DALIGNER_OPT} -r ${COR_CONTIG_DALIGNER_RUNID}"
    fi
    if [[ -n ${COR_CONTIG_DALIGNER_OLEN} ]]
    then
        CONTIG_DALIGNER_OPT="${CONTIG_DALIGNER_OPT} -l ${COR_CONTIG_DALIGNER_OLEN}"
    fi    
    if [[ -n ${COR_CONTIG_DALIGNER_MEM} && ${COR_CONTIG_DALIGNER_MEM} -gt 0 ]]
    then
        CONTIG_DALIGNER_OPT="${CONTIG_DALIGNER_OPT} -M ${COR_CONTIG_DALIGNER_MEM}"
    fi    
    if [[ -n ${COR_CONTIG_DALIGNER_HITS} ]]
    then
        CONTIG_DALIGNER_OPT="${CONTIG_DALIGNER_OPT} -h ${COR_CONTIG_DALIGNER_HITS}"
    fi        
    if [[ -n ${COR_CONTIG_DALIGNER_T} ]]
    then
        CONTIG_DALIGNER_OPT="${CONTIG_DALIGNER_OPT} -t ${COR_CONTIG_DALIGNER_T}"
    fi  
    if [[ -n ${COR_CONTIG_DALIGNER_MASK} ]]
    then
        for x in ${COR_CONTIG_DALIGNER_MASK}
        do 
            CONTIG_DALIGNER_OPT="${CONTIG_DALIGNER_OPT} -m ${x}"
        done
    fi
    if [[ -n ${COR_CONTIG_DALIGNER_TRACESPACE} && ${COR_CONTIG_DALIGNER_TRACESPACE} -gt 0 ]]
    then
        CONTIG_DALIGNER_OPT="${CONTIG_DALIGNER_OPT} -s ${COR_CONTIG_DALIGNER_TRACESPACE}"
    fi
    if [[ -n ${THREADS_daligner} ]]
    then 
        CONTIG_DALIGNER_OPT="${CONTIG_DALIGNER_OPT} -j ${THREADS_daligner}"
    fi
}

function setLAmergeOptions()
{
    CONTIG_LAMERGE_OPT=""
    if [[ -n ${COR_CONTIG_LAMERGE_NFILES} && ${COR_CONTIG_LAMERGE_NFILES} -gt 0 ]]
    then
        CONTIG_LAMERGE_OPT="${CONTIG_LAMERGE_OPT} -n ${COR_CONTIG_LAMERGE_NFILES}"
    fi
}    

function setLArepeatOptions()
{
	if [[ ${#COR_CONTIG_LAREPEAT_LEAVE_COV[*]} -ne ${#COR_CONTIG_LAREPEAT_ENTER_COV[*]} || ${#COR_CONTIG_LAREPEAT_ENTER_COV[*]} -ne ${#COR_CONTIG_LAREPEAT_COV[*]} ]]
    then 
        (>&2 echo "LArepeat number of elements of COR_CONTIG_LAREPEAT_LEAVE_COV and COR_CONTIG_LAREPEAT_ENTER_COV and COR_CONTIG_LAREPEAT_COV differ")
        (>&2 echo "they must be of the same length")
        exit 1
    fi

    numRepeatTracks=${#COR_CONTIG_LAREPEAT_LEAVE_COV[*]}

    # define array variable - because we may want to create several repeat tracks in one run
    unset CONTIG_LAREPEAT_OPT
    ### find and set LArepeat options     
    
    for x in $(seq 0 $((${numRepeatTracks}-1)))
    do 
        tmp=""
        tmp="${tmp} -l ${COR_CONTIG_LAREPEAT_LEAVE_COV[$x]}"
        tmp="${tmp} -h ${COR_CONTIG_LAREPEAT_ENTER_COV[$x]}"

        if [[ -n ${COR_CONTIG_LAREPEAT_OLEN} && ${COR_CONTIG_LAREPEAT_OLEN} -gt 0 ]]
        then
            tmp="${tmp} -o ${COR_CONTIG_LAREPEAT_OLEN}"
        fi

        if [[ -n ${COR_CONTIG_LAREPEAT_IDENTITYOVLS} && ${COR_CONTIG_LAREPEAT_IDENTITYOVLS} -gt 0 ]]
        then
            tmp="${tmp} -I"
        fi                

        tmp="${tmp} -c ${COR_CONTIG_LAREPEAT_COV[$x]}"
        tmp="${tmp} -t repeats_c${COR_CONTIG_LAREPEAT_COV[$x]}_l${COR_CONTIG_LAREPEAT_LEAVE_COV[$x]}h${COR_CONTIG_LAREPEAT_ENTER_COV[$x]}"
        CONTIG_LAREPEAT_OPT[$x]=${tmp}
    done 
}

function setTKmergeOptions() 
{
    CONTIG_TKMERGE_OPT=""
    if [[ -n ${COR_CONTIG_TKMERGE_DELETE} && ${COR_CONTIG_TKMERGE_DELETE} -ne 0 ]]
    then
        CONTIG_TKMERGE_OPT="${CONTIG_TKMERGE_OPT} -d"
    fi
}

function setTKcombineOptions() 
{
    ignoreDelete=$1
    CONTIG_TKCOMBINE_OPT=""
    if [[ ${ignoreDelete} -eq 0 && -n ${COR_CONTIG_TKCOMBINE_DELETE} && ${COR_CONTIG_TKCOMBINE_DELETE} -ne 0 ]]
    then
        CONTIG_TKCOMBINE_OPT="${CONTIG_TKCOMBINE_OPT} -d"
    fi
    if [[ -n ${COR_CONTIG_TKCOMBINE_VERBOSE} && ${COR_CONTIG_TKCOMBINE_VERBOSE} -ne 0 ]]
    then
        CONTIG_TKCOMBINE_OPT="${CONTIG_TKCOMBINE_OPT} -v"
    fi
}

function setCTanalyzeOptions()
{
	CONTIG_CTANALYZE_OPT=""
	if [[ -n ${COR_CONTIG_CTANALYZE_VERBOSE} && ${COR_CONTIG_CTANALYZE_VERBOSE} -ne 0 ]]
    then
        CONTIG_CTANALYZE_OPT="${CONTIG_CTANALYZE_OPT} -v"
    fi
        
    if [[ -n ${COR_CONTIG_CTANALYZE_EXPRAWREADCOV} && ${COR_CONTIG_CTANALYZE_EXPRAWREADCOV} -gt 0 ]]
    then
        CONTIG_CTANALYZE_OPT="${CONTIG_CTANALYZE_OPT} -c ${COR_CONTIG_CTANALYZE_EXPRAWREADCOV}"
    fi
    
    if [[ -n ${COR_CONTIG_CTANALYZE_MAXSPURLEN} && ${COR_CONTIG_CTANALYZE_MAXSPURLEN} -gt 0 ]]
    then
        CONTIG_CTANALYZE_OPT="${CONTIG_CTANALYZE_OPT} -s ${COR_CONTIG_CTANALYZE_MAXSPURLEN}"
    fi
    
    if [[ -n ${COR_CONTIG_CTANALYZE_MAXTIPLEN} && ${COR_CONTIG_CTANALYZE_MAXTIPLEN} -gt 0 ]]
    then
        CONTIG_CTANALYZE_OPT="${CONTIG_CTANALYZE_OPT} -e ${COR_CONTIG_CTANALYZE_MAXTIPLEN}"
    fi	
    	    
    if [[ -n ${COR_CONTIG_CTANALYZE_TRIMTRACK_PATCHEDREADS} ]]
    then
        CONTIG_CTANALYZE_OPT="${CONTIG_CTANALYZE_OPT} -t ${COR_CONTIG_CTANALYZE_TRIMTRACK_PATCHEDREADS}"
    fi
    
    if [[ -n ${COR_CONTIG_CTANALYZE_DIR} ]]
    then
    	if [[ -d ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${COR_CONTIG_CTANALYZE_DIR} ]]; 
    	then 
    		mv ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${COR_CONTIG_CTANALYZE_DIR} ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${COR_CONTIG_CTANALYZE_DIR}_$(date '+%Y-%m-%d_%H-%M-%S'); 
    	fi 
    	mkdir -p ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${COR_CONTIG_CTANALYZE_DIR}
        CONTIG_CTANALYZE_OPT="${CONTIG_CTANALYZE_OPT} -d ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${COR_CONTIG_CTANALYZE_DIR}"
    fi
    
    if [[ -n ${COR_CONTIG_CTANALYZE_CONTIGREPEATTRACK} ]]
    then
        CONTIG_CTANALYZE_OPT="${CONTIG_CTANALYZE_OPT} -r ${COR_CONTIG_CTANALYZE_CONTIGREPEATTRACK}"
    fi
    
    if [[ -n ${COR_CONTIG_CTANALYZE_READREPEATTRACK} ]]
    then
        CONTIG_CTANALYZE_OPT="${CONTIG_CTANALYZE_OPT} -R ${COR_CONTIG_CTANALYZE_READREPEATTRACK}"
    fi	

    if [[ -n ${COR_CONTIG_CTANALYZE_FUZZYSVLEN} ]]
    then
        CONTIG_CTANALYZE_OPT="${CONTIG_CTANALYZE_OPT} -f ${COR_CONTIG_CTANALYZE_FUZZYSVLEN}"
    fi	
    
    if [[ -n ${COR_CONTIG_CTANALYZE_MINPRIMLEN} ]]
    then
        CONTIG_CTANALYZE_OPT="${CONTIG_CTANALYZE_OPT} -L ${COR_CONTIG_CTANALYZE_MINPRIMLEN}"
    fi
    
    if [[ -n ${COR_CONTIG_CTANALYZE_MINPRIMCREADS} ]]
    then
        CONTIG_CTANALYZE_OPT="${CONTIG_CTANALYZE_OPT} -N ${COR_CONTIG_CTANALYZE_MINPRIMCREADS}"
    fi
    
    if [[ -n ${COR_CONTIG_CTANALYZE_MAXREPEATPERC} ]]
    then
        CONTIG_CTANALYZE_OPT="${CONTIG_CTANALYZE_OPT} -P ${COR_CONTIG_CTANALYZE_MAXREPEATPERC}"
    fi
}

function setDalignerOptions()
{
    CONTIG_DALIGNER_OPT=""
    if [[ -n ${COR_CONTIG_DALIGNER_IDENTITY_OVLS} && ${COR_CONTIG_DALIGNER_IDENTITY_OVLS} -gt 0 ]]
    then
        CONTIG_DALIGNER_OPT="${CONTIG_DALIGNER_OPT} -I"
    fi
    if [[ -n ${COR_CONTIG_DALIGNER_KMER} && ${COR_CONTIG_DALIGNER_KMER} -gt 0 ]]
    then
        CONTIG_DALIGNER_OPT="${CONTIG_DALIGNER_OPT} -k ${COR_CONTIG_DALIGNER_KMER}"
    fi
    if [[ -n ${COR_CONTIG_DALIGNER_ERR} ]]
    then
        CONTIG_DALIGNER_OPT="${CONTIG_DALIGNER_OPT} -e ${COR_CONTIG_DALIGNER_ERR}"
    fi
    if [[ -n ${COR_CONTIG_DALIGNER_BIAS} && ${COR_CONTIG_DALIGNER_BIAS} -eq 1 ]]
    then
        CONTIG_DALIGNER_OPT="${CONTIG_DALIGNER_OPT} -b"
    fi
    if [[ -n ${COR_CONTIG_DALIGNER_VERBOSE} && ${COR_CONTIG_DALIGNER_VERBOSE} -ne 0 ]]
    then
        CONTIG_DALIGNER_OPT="${CONTIG_DALIGNER_OPT} -v"
    fi
    if [[ -n ${COR_CONTIG_DALIGNER_TRACESPACE} && ${COR_CONTIG_DALIGNER_TRACESPACE} -gt 0 ]]
    then
        CONTIG_DALIGNER_OPT="${CONTIG_DALIGNER_OPT} -s ${COR_CONTIG_DALIGNER_TRACESPACE}"
    fi
    if [[ -n ${COR_CONTIG_DALIGNER_RUNID} ]]
    then
        CONTIG_DALIGNER_OPT="${CONTIG_DALIGNER_OPT} -r ${COR_CONTIG_DALIGNER_RUNID}"
    fi
    if [[ -n ${COR_CONTIG_DALIGNER_T} ]]
    then
        CONTIG_DALIGNER_OPT="${CONTIG_DALIGNER_OPT} -t ${COR_CONTIG_DALIGNER_T}"
    fi  
    if [[ -n ${COR_CONTIG_DALIGNER_ASYMMETRIC} && ${COR_CONTIG_DALIGNER_ASYMMETRIC} -ne 0 ]]
    then
        CONTIG_DALIGNER_OPT="${CONTIG_DALIGNER_OPT} -A"
    fi    
    if [[ -n ${COR_CONTIG_DALIGNER_MEM} && ${COR_CONTIG_DALIGNER_MEM} -ne 0 ]]
    then
        CONTIG_DALIGNER_OPT="${CONTIG_DALIGNER_OPT} -M ${COR_CONTIG_DALIGNER_MEM}"
    fi    
    if [[ -n ${COR_CONTIG_DALIGNER_MASK} ]]
    then
        for x in ${COR_CONTIG_DALIGNER_MASK}
        do 
            CONTIG_DALIGNER_OPT="${CONTIG_DALIGNER_OPT} -m ${x}"
        done
    fi
    if [[ -n ${THREADS_daligner} ]]
    then 
        CONTIG_DALIGNER_OPT="${CONTIG_DALIGNER_OPT} -j ${THREADS_daligner}"
    fi
    if [ ! -n ${COR_CONTIG_DALIGNER_DAL} ]
    then
        COR_CONTIG_DALIGNER_DAL=8
    fi 
}

if [[ -z ${COR_DIR} ]]
then 
    COR_DIR=correction
fi

if [[ -z ${ANALYZE_DIR} ]]
then 
    ANALYZE_DIR=analyze
fi

### set filter directory

if [[ -z ${FIX_FILT_OUTDIR} ]]
then
    FIX_FILT_OUTDIR="m1"
fi

if [[ -f ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${CONT_DB%.db}.db ]]
then
	contigblocks=$(getNumOfDbBlocks ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${CONT_DB%.db}.db)
fi

## ensure some paths
if [[ -z "${MARVEL_SOURCE_PATH}" ]]
then 
    (>&2 echo "ERROR - You have to set MARVEL_SOURCE_PATH. Used to report git version.")
    exit 1
fi

if [[ -z "${DAZZLER_SOURCE_PATH}" ]]
then 
    (>&2 echo "ERROR - You have to set DAZZLER_SOURCE_PATH. Used to report git version.")
    exit 1
fi

if [[ -z "${DACCORD_SOURCE_PATH}" ]]
then 
    (>&2 echo "ERROR - You have to set DACCORD_SOURCE_PATH. Used to report git version.")
    exit 1
fi

myTypes=("1-createCorrectedContigDB, 2-DBdust, 3-Catrack, 4-datander, 5-TANmask, 6-Catrack, 7-daligner, 08_LAmerge, 09_LArepeat, 10_TKmerge, 11_TKcombine, 12_LAfilter, 13_LAmerge, 14_CTanalyze, 15_CTstatistics")
#type-0 steps: 01_createCorrectedContigDB, 02_DBdust, 03_Catrack, 04_datander, 05_TANmask, 06_Catrack, 07_daligner, 08_LAmerge, 09_LArepeat, 10_TKmerge, 11_TKcombine, 12_LAfilter, 13_LAmerge, 14_CTanalyze, 15_CTstatistics
if [[ ${COR_CONTIG_TYPE} -eq 0 ]]
then 
    ### createCorrectedContigDB
    if [[ ${currentStep} -eq 1 ]]
    then
        ### clean up plans 
        for x in $(ls cont_01_*_*_${CONT_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        setcreateCorrectedContigDBOptions
        
        # check created fasta files
        if [[ ! -d ${FIX_FILT_OUTDIR}/${COR_DIR}/contigs ]]
        then
        	(>&2 echo "missing directory ${FIX_FILT_OUTDIR}/${COR_DIR}/contigs")
    		exit 1	
        fi
        
        echo "if [[ -d ${FIX_FILT_OUTDIR}/${ANALYZE_DIR} ]]; then mv ${FIX_FILT_OUTDIR}/${ANALYZE_DIR} ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}_$(date '+%Y-%m-%d_%H-%M-%S'); fi && mkdir -p ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}" > cont_01_prepDB_single_${CONT_DB%.db}.${slurmID}.plan
        # link FIX- RAW- and COR-DB into working directory 
        echo "ln -s -r .${FIX_DB%db}.* ${FIX_DB%db}.db ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}" >> cont_01_prepDB_single_${CONT_DB%.db}.${slurmID}.plan
        echo "ln -s -r ${FIX_FILT_OUTDIR}/${COR_DIR}/.${COR_DB%db}.* ${FIX_FILT_OUTDIR}/${COR_DIR}/${COR_DB%db}.db ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}" >> cont_01_prepDB_single_${CONT_DB%.db}.${slurmID}.plan
		echo "ln -s -r .${RAW_DB%db}.* ${RAW_DB%db}.db ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}" >> cont_01_prepDB_single_${CONT_DB%.db}.${slurmID}.plan
		
		first=1
		for x in ${FIX_FILT_OUTDIR}/${COR_DIR}/contigs/*.fasta
		do 
			if [[ ! -f ${x} ]]
			then 
				(>&2 echo "missing corrected contigs $x in directory ${FIX_FILT_OUTDIR}/${COR_DIR}/contigs")
    			exit 1
			fi	
			if [[  $first -eq 1 ]]
			then
				first=0 
				echo "${MARVEL_PATH}/bin/FA2db -x0 -c path -c ends -c length -c rreads -c preads -c creads ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${CONT_DB} ${x}"
				echo "${MARVEL_PATH}/bin/DBsplit ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${CONT_DB}"
			else
				echo "${MARVEL_PATH}/bin/FA2db -x0 -c path -c ends -c length -c rreads -c preads -c creads ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${CONT_DB} ${x}"				
			fi				
		done >> cont_01_prepDB_single_${CONT_DB%.db}.${slurmID}.plan
		
		echo "mkdir -p ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/correctedContigs_dazzler" >> cont_01_prepDB_single_${CONT_DB%.db}.${slurmID}.plan
		
		## create a proper dazzler fasta header
        for x in ${FIX_FILT_OUTDIR}/${COR_DIR}/contigs/*.fasta
        do  
            echo "${DACCORD_PATH}/bin/fastaidrename < ${x} | awk '{print \$1}' > ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/correctedContigs_dazzler/$(basename ${x%.fasta})_dazzler.fasta"            
		done >> cont_01_prepDB_single_${CONT_DB%.db}.${slurmID}.plan        

        # create dazzler db
    	first=1 
        for x in ${FIX_FILT_OUTDIR}/${COR_DIR}/contigs/*.fasta
        do
            if [[ ${first} -eq 1 ]]
	        then
	        	first=0
				echo "${DAZZLER_PATH}/bin/fasta2DB -v ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${CONT_DAZZ_DB} ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/correctedContigs_dazzler/$(basename ${x%.fasta})_dazzler.fasta"
				echo "${DAZZLER_PATH}/bin/DBsplit ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${CONT_DAZZ_DB}"	        		
	        else
	        	echo "${DAZZLER_PATH}/bin/fasta2DB -v ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${CONT_DAZZ_DB} ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/correctedContigs_dazzler/$(basename ${x%.fasta})_dazzler.fasta"
	        fi              
		done >> cont_01_prepDB_single_${CONT_DB%.db}.${slurmID}.plan
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > cont_01_prepDB_single_${CONT_DB%.db}.${slurmID}.version
		echo "DAZZ_DB $(git --git-dir=${DAZZLER_SOURCE_PATH}/DAZZ_DB/.git rev-parse --short HEAD)" >> cont_01_prepDB_single_${CONT_DB%.db}.${slurmID}.version       
    ### DBdust
    elif [[ ${currentStep} -eq 2 ]]
    then
        ### clean up plans 
        for x in $(ls cont_02_*_*_${CONT_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        ### find and set DBdust options 
        setDBdustOptions
        ### create DBdust commands 
        for x in $(seq 1 ${contigblocks})
        do 
            echo "${MARVEL_PATH}/bin/DBdust${CONTIG_DBDUST_OPT} ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${CONT_DB%.db}.${x}"
            echo "${DAZZLER_PATH}/bin/DBdust${CONTIG_DBDUST_OPT} ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${CONT_DAZZ_DB%.db}.${x}"
		done > cont_02_DBdust_block_${CONT_DB%.db}.${slurmID}.plan
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > cont_02_DBdust_block_${CONT_DB%.db}.${slurmID}.version
		echo "DAZZ_DB $(git --git-dir=${DAZZLER_SOURCE_PATH}/DAZZ_DB/.git rev-parse --short HEAD)" >> cont_02_DBdust_block_${CONT_DB%.db}.${slurmID}.version
	##$ Catrack
 	elif [[ ${currentStep} -eq 3 ]]
    then 
        ### clean up plans 
        for x in $(ls cont_03_*_*_${CONT_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        
        ### find and set Catrack options 
        setCatrackOptions
        ### create Catrack command
        echo "${MARVEL_PATH}/bin/Catrack${CONTIG_CATRACK_OPT} ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${CONT_DB%.db} dust" > cont_03_Catrack_single_${CONT_DB%.db}.${slurmID}.plan
        echo "${DAZZLER_PATH}/bin/Catrack${CONTIG_CATRACK_OPT} ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${CONT_DAZZ_DB%.db} dust" >> cont_03_Catrack_single_${CONT_DB%.db}.${slurmID}.plan
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > cont_03_Catrack_single_${CONT_DB%.db}.${slurmID}.version
		echo "DAZZ_DB $(git --git-dir=${DAZZLER_SOURCE_PATH}/DAZZ_DB/.git rev-parse --short HEAD)" >> cont_03_Catrack_single_${CONT_DB%.db}.${slurmID}.version       
    ### datander
    elif [[ ${currentStep} -eq 4 ]]
    then 
        ### clean up plans 
        for x in $(ls cont_04_*_*_${CONT_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done     
        ### find and set datander options 
        setDatanderOptions
        d=$(pwd)
        ### create datander commands
        for x in $(seq 1 ${contigblocks})
        do 
            echo "cd ${FIX_FILT_OUTDIR}/${ANALYZE_DIR} && ${MARVEL_PATH}/bin/datander${CONTIG_DATANDER_OPT} ${CONT_DB%.db}.${x} && cd ${d}"
		done > cont_04_datander_block_${CONT_DB%.db}.${slurmID}.plan
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > cont_04_datander_block_${CONT_DB%.db}.${slurmID}.version
	### TANmask
    elif [[ ${currentStep} -eq 5 ]]
    then 
        ### clean up plans 
        for x in $(ls cont_05_*_*_${CONT_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done     
        ### find and set TANmask options 
        setTANmaskOptions
        ### create TANmask commands
        for x in $(seq 1 ${contigblocks})
        do 
            echo "${MARVEL_PATH}/bin/TANmask${CONTIG_TANMASK_OPT} ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${CONT_DB%.db} ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${COR_CONTIG_DATANDER_FOLDER}/${CONT_DB%.db}.${x}.${CONT_DB%.db}.${x}.las" 
		done > cont_05_TANmask_block_${CONT_DB%.db}.${slurmID}.plan
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > cont_05_TANmask_block_${CONT_DB%.db}.${slurmID}.version
	### Catrack	
    elif [[ ${currentStep} -eq 6 ]]
    then 
        ### clean up plans 
        for x in $(ls cont_06_*_*_${CONT_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        ### find and set Catrack options
        if [[ -z ${CONTIG_CATRACK_OPT} ]] 
        then
            setCatrackOptions
        fi
        ### create Catrack command
        echo "${MARVEL_PATH}/bin/Catrack${CONTIG_CATRACK_OPT} ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${CONT_DB%.db} ${COR_CONTIG_TANMASK_TRACK}" > cont_06_Catrack_single_${CONT_DB%.db}.${slurmID}.plan
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > cont_06_Catrack_single_${CONT_DB%.db}.${slurmID}.version
   	### daligner 
    elif [[ ${currentStep} -eq 7 ]]
    then
        ### clean up plans 
        for x in $(ls cont_07_*_*_${CONT_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        ### find and set daligner options 
        setDalignerOptions
        cmdLine=1
        d=$(pwd)
        ### create daligner commands
        for x in $(seq 1 ${contigblocks})
        do 
            if [[ -n ${COR_CONTIG_DALIGNER_NUMACTL} && ${COR_CONTIG_DALIGNER_NUMACTL} -gt 0 ]] && [[ "x${SLURM_NUMACTL}" == "x" || ${SLURM_NUMACTL} -eq 0 ]]
            then
                if [[ $((${cmdLine} % 2)) -eq  0 ]]
                then
                    NUMACTL="numactl -m0 -N0 "
                else
                    NUMACTL="numactl -m1 -N1 "    
                fi
            else
                NUMACTL=""
            fi
            cmd="cd ${FIX_FILT_OUTDIR}/${ANALYZE_DIR} && ${NUMACTL}${MARVEL_PATH}/bin/daligner${CONTIG_DALIGNER_OPT} ${CONT_DB%.db}.${x}"
            cmdLine=$((${cmdLine}+1))
            count=0
            for y in $(seq ${x} ${contigblocks})
            do  
            	if [[ $count -lt ${COR_CONTIG_DALIGNER_DAL} || ${count} -lt 4 ]]
                then
                    cmd="${cmd} ${CONT_DB%.db}.${y}"
                    count=$((${count}+1))
                else    
                    echo "${cmd} && cd $d"
                    if [[ -n ${COR_CONTIG_DALIGNER_NUMACTL} && ${COR_CONTIG_DALIGNER_NUMACTL} -gt 0 ]] && [[ "x${SLURM_NUMACTL}" == "x" || ${SLURM_NUMACTL} -eq 0 ]]
                    then
                        if [[ $((${cmdLine} % 2)) -eq  0 ]]
                        then
                            NUMACTL="numactl -m0 -N0 "
                        else
                            NUMACTL="numactl -m1 -N1 "    
                        fi
                    else
                        NUMACTL=""
                    fi
                    cmd="cd ${FIX_FILT_OUTDIR}/${ANALYZE_DIR} && ${NUMACTL}${MARVEL_PATH}/bin/daligner${CONTIG_DALIGNER_OPT} ${CONT_DB%.db}.${x} ${CONT_DB%.db}.${y}"
                    cmdLine=$((${cmdLine}+1))
                    count=1
                fi
            done
            echo "${cmd} && cd ${d}"
		done > cont_07_daligner_block_${CONT_DB%.db}.${slurmID}.plan
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > cont_07_daligner_block_${CONT_DB%.db}.${slurmID}.version
    #### LAmerge
    elif [[ ${currentStep} -eq 8 ]]
    then
        ### clean up plans 
        for x in $(ls cont_08_*_*_${CONT_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        
        ### find and set LAmerge options 
        setLAmergeOptions
        ### create LAmerge commands
        for x in $(seq 1 ${contigblocks}); 
        do  
            echo "${MARVEL_PATH}/bin/LAmerge${CONTIG_LAMERGE_OPT} ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${CONT_DB%.db} ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${CONT_DB%.db}.${x}.las ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/$(getSubDirName ${COR_CONTIG_DALIGNER_RUNID} ${x})"
		done > cont_08_LAmerge_block_${CONT_DB%.db}.${slurmID}.plan
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > cont_08_LAmerge_block_${CONT_DB%.db}.${slurmID}.version   
    #### LArepeat
    elif [[ ${currentStep} -eq 9 ]]
    then
        ### clean up plans 
        for x in $(ls cont_09_*_*_${CONT_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 

        setLArepeatOptions
        
        if [[ ${numRepeatTracks} -eq 0 ]]
        then 
            exit 1
        fi    
    
        for x in $(seq 0 $((${numRepeatTracks}-1)))
        do 
            ### create LArepeat commands
            for y in $(seq 1 ${contigblocks})
            do 
                echo "${MARVEL_PATH}/bin/LArepeat${CONTIG_LAREPEAT_OPT[$x]} -b ${y} ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${CONT_DB%.db} ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${CONT_DB%.db}.${y}.las"
            done
		done > cont_09_LArepeat_block_${CONT_DB%.db}.${slurmID}.plan
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > cont_09_LArepeat_block_${CONT_DB%.db}.${slurmID}.version      
    #### TKmerge 
    elif [[ ${currentStep} -eq 10 ]]
    then
        ### clean up plans 
        for x in $(ls cont_10_*_*_${CONT_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done         
        # we need the name of the repeat track, especially if the plan starts with step10
        if [[ -z ${CONTIG_LAREPEAT_OPT} ]]
        then 
            setLArepeatOptions
        fi
        ### find and set TKmerge options 
        if [[ -z ${CONTIG_TKMERGE_OPT} ]]
        then 
            setTKmergeOptions
        fi
        
        if [[ ${numRepeatTracks} -eq 0 ]]
        then 
            exit 1
        fi 
        ### create TKmerge command
        for x in $(seq 0 $((${numRepeatTracks}-1)))
        do 
            echo "${MARVEL_PATH}/bin/TKmerge${CONTIG_TKMERGE_OPT} ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${CONT_DB%.db} $(echo ${CONTIG_LAREPEAT_OPT[${x}]} | awk '{print $NF}')" 
		done > cont_10_TKmerge_single_${CONT_DB%.db}.${slurmID}.plan  
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > cont_10_TKmerge_single_${CONT_DB%.db}.${slurmID}.version     
    #### TKcombine 
    elif [[ ${currentStep} -eq 11 ]]
    then
        ### clean up plans 
        for x in $(ls cont_11_*_*_${CONT_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done             
         # we need the name of the repeat track, especially if the plan starts with step11
        if [[ -z ${CONTIG_LAREPEAT_OPT} ]]
        then 
            setLArepeatOptions
        fi
        ### find and set TKcombine options
        setTKcombineOptions 1

        if [[ ${numRepeatTracks} -eq 0 ]]
        then 
            exit 1
        fi 
        ### create TKmerge command
        for x in $(seq 0 $((${numRepeatTracks}-1)))
        do 
            tmp=$(echo ${CONTIG_LAREPEAT_OPT[${x}]} | awk '{print $NF}')
            echo "${MARVEL_PATH}/bin/TKcombine${CONTIG_TKCOMBINE_OPT} ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${CONT_DB%.db} ${tmp}_${COR_CONTIG_TANMASK_TRACK} ${tmp} ${COR_CONTIG_TANMASK_TRACK}"
            echo "${MARVEL_PATH}/bin/TKcombine${CONTIG_TKCOMBINE_OPT} ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${CONT_DB%.db} ${tmp}_${COR_CONTIG_TANMASK_TRACK}_dust ${tmp}_${COR_CONTIG_TANMASK_TRACK} dust" 
		done > cont_11_TKcombine_single_${CONT_DB%.db}.${slurmID}.plan
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > cont_11_TKcombine_single_${CONT_DB%.db}.${slurmID}.version         
	### LAfilterChains
    elif [[ ${currentStep} -eq 12 ]]
    then
        ### clean up plans 
        for x in $(ls cont_12_*_*_${CONT_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
 
        if [[ -z ${CONTIG_LAFILTERCHAINS_OPT} ]]
        then 
            setLAfilterChainsOptions
        fi

        ### create LAfilter commands
        contigblocks=$(getNumOfDbBlocks ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${CONT_DB%.db}.db)       
        for x in $(seq 1 ${contigblocks});
        do 
            echo "${MARVEL_PATH}/bin/LAfilterChains${CONTIG_LAFILTERCHAINS_OPT} ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${CONT_DB%.db} ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${CONT_DB%.db}.${x}.las ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${CONT_DB%.db}.${x}.filt.las"
		done > cont_12_LAfilter_block_${CONT_DB%.db}.${slurmID}.plan 
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > cont_12_LAfilter_block_${CONT_DB%.db}.${slurmID}.version
	### LAmerge
    elif [[ ${currentStep} -eq 13 ]]
    then
        ### clean up plans 
        for x in $(ls cont_13_*_*_${CONT_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        
        ### find and set LAmerge options 
        setLAmergeOptions
        
        echo "${MARVEL_PATH}/bin/LAmerge${CONTIG_LAMERGE_OPT} -S filt ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${CONT_DB%.db} ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${CONT_DB%.db}.filt.las" > cont_13_LAmerge_single_${CONT_DB%.db}.${slurmID}.plan
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > cont_13_LAmerge_single_${CONT_DB%.db}.${slurmID}.version
	### CTanalyze
	elif [[ ${currentStep} -eq 14 ]]
    then
		### clean up plans 
        for x in $(ls cont_14_*_*_${CONT_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        
        ### find and set CTanalyze options 
        setCTanalyzeOptions
		echo "${MARVEL_PATH}/bin/CTanalyze${CONTIG_CTANALYZE_OPT} -C ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${CONT_DB%.db} ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${CONT_DB%.db}.filt.las -F ${FIX_FILT_OUTDIR}/${FIX_DB%.DB} ${FIX_FILT_OUTDIR}/${FIX_DB%.DB}.filt.las -D ${FIX_FILT_OUTDIR}/${COR_DIR}/${COR_DB%.db}" > cont_14_CTanalyze_single_${CONT_DB%.db}.${slurmID}.plan
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > cont_14_CTanalyze_single_${CONT_DB%.db}.${slurmID}.version
    ### 15_CTstatistics 
	elif [[ ${currentStep} -eq 15 ]]
    then
		### clean up plans 
        for x in $(ls cont_15_*_*_${CONT_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
                
        ### run slurm stats - on the master node !!! Because sacct is not available on compute nodes
        if [[ -n ${COR_CONTIG_CTANALYZE_FULLSTATS} && ${COR_CONTIG_CTANALYZE_FULLSTATS} -gt 0 ]]
        then
	        if [[ $(hostname) == "falcon1" || $(hostname) == "falcon2" ]]
	        then 
	        	bash ${SUBMIT_SCRIPTS_PATH}/slurmStats.sh ${configFile}
	    	else
	        	cwd=$(pwd)
	        	ssh falcon "cd ${cwd} && bash ${SUBMIT_SCRIPTS_PATH}/slurmStats.sh ${configFile}"
	    	fi
		fi
        ### create assemblyStats plan 
        echo "${SUBMIT_SCRIPTS_PATH}/assemblyStats.sh ${configFile} 8" > cont_15_CTstatistics_single_${CONT_DB%.db}.${slurmID}.plan
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > cont_15_CTstatistics_single_${CONT_DB%.db}.${slurmID}.version 
	else
        (>&2 echo "step ${currentStep} in COR_CONTIG_TYPE ${COR_CONTIG_TYPE} not supported")
        (>&2 echo "valid steps are: ${myTypes[${RAW_REPMASK_TYPE}]}")
        exit 1            
    fi
else
    (>&2 echo "unknown COR_CONTIG_TYPE ${COR_CONTIG_TYPE}")
    (>&2 echo "supported types")
    x=0; while [ $x -lt ${#myTypes[*]} ]; do (>&2 echo "${myTypes[${x}]}"); done 
    
    exit 1
fi

exit 0