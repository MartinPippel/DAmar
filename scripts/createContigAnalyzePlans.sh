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

function setLAfilterOptions()
{
    CONTIG_LAFILTER_OPT=""

    if [[ -n ${COR_CONTIG_LAFILTER_NREP} && ${COR_CONTIG_LAFILTER_NREP} -ne 0 ]]
    then
        CONTIG_LAFILTER_OPT="${CONTIG_LAFILTER_OPT} -n ${COR_CONTIG_LAFILTER_NREP}"
    fi
    if [[ -n ${COR_CONTIG_LAFILTER_VERBOSE} && ${COR_CONTIG_LAFILTER_VERBOSE} -ne 0 ]]
    then
        CONTIG_LAFILTER_OPT="${CONTIG_LAFILTER_OPT} -v"
    fi
    if [[ -n ${COR_CONTIG_LAFILTER_PURGE} && ${COR_CONTIG_LAFILTER_PURGE} -ne 0 ]]
    then
        CONTIG_LAFILTER_OPT="${CONTIG_LAFILTER_OPT} -p"
    fi
    if [[ -n ${COR_CONTIG_LAFILTER_OLEN} && ${COR_CONTIG_LAFILTER_OLEN} -ne 0 ]]
    then
        CONTIG_LAFILTER_OPT="${CONTIG_LAFILTER_OPT} -o ${COR_CONTIG_LAFILTER_OLEN}"
    fi    
    if [[ -n ${COR_CONTIG_LAFILTER_RLEN} && ${COR_CONTIG_LAFILTER_RLEN} -ne 0 ]]
    then
        CONTIG_LAFILTER_OPT="${CONTIG_LAFILTER_OPT} -l ${COR_CONTIG_LAFILTER_RLEN}"
    fi   

    if [[ -n ${COR_CONTIG_LAFILTER_DIF} && ${COR_CONTIG_LAFILTER_DIF} -ne 0 ]]
    then
        CONTIG_LAFILTER_OPT="${CONTIG_LAFILTER_OPT} -d ${COR_CONTIG_LAFILTER_DIF}"
    fi

    if [[ -n ${COR_CONTIG_LAFILTER_UBAS} ]]
    then
        CONTIG_LAFILTER_OPT="${CONTIG_LAFILTER_OPT} -u ${COR_CONTIG_LAFILTER_UBAS}"
    fi
    if [[ -n ${COR_CONTIG_LAFILTER_PRELOAD} && ${COR_CONTIG_LAFILTER_PRELOAD} -ne 0 ]]
    then
        CONTIG_LAFILTER_OPT="${CONTIG_LAFILTER_OPT} -L"
    fi    
    if [[ -n ${COR_CONTIG_LAFILTER_MERGEREPEATS} && ${COR_CONTIG_LAFILTER_MERGEREPEATS} -ne 0 ]]
    then
        CONTIG_LAFILTER_OPT="${CONTIG_LAFILTER_OPT} -y ${COR_CONTIG_LAFILTER_MERGEREPEATS}"
    fi    
    if [[ -n ${COR_CONTIG_LAFILTER_MERGEREPEATTIPS} && ${COR_CONTIG_LAFILTER_MERGEREPEATTIPS} -ne 0 ]]
    then
        CONTIG_LAFILTER_OPT="${CONTIG_LAFILTER_OPT} -Y ${COR_CONTIG_LAFILTER_MERGEREPEATTIPS}"
    fi    
    if [[ -n ${COR_CONTIG_LAFILTER_MINTIPCOV} && ${COR_CONTIG_LAFILTER_MINTIPCOV} -gt 0 ]]
    then
        CONTIG_LAFILTER_OPT="${CONTIG_LAFILTER_OPT} -z ${COR_CONTIG_LAFILTER_MINTIPCOV}"
    fi            
    if [[ -n ${COR_CONTIG_LAFILTER_MULTIMAPPER} && ${COR_CONTIG_LAFILTER_MULTIMAPPER} -gt 0 ]]
    then
        if [[ ${COR_CONTIG_LAFILTER_MULTIMAPPER} -eq 1 ]]
        then
            CONTIG_LAFILTER_OPT="${CONTIG_LAFILTER_OPT} -w"
        else
            CONTIG_LAFILTER_OPT="${CONTIG_LAFILTER_OPT} -W"
        fi
    fi            

    if [[ -n ${COR_CONTIG_LAFILTER_RESOLVE_REPEATS} && ${COR_CONTIG_LAFILTER_RESOLVE_REPEATS} -gt 0 && ${COR_CONTIG_LAFILTER_RESOLVE_REPEATS} -lt 4 ]]
    then
        tmp=""
        mode="m"
        if [[ -n ${COR_CONTIG_LAFILTER_RESOLVE_REPEATS_AGG} && ${COR_CONTIG_LAFILTER_RESOLVE_REPEATS_AGG} -ne 0 ]]
        then
            mode="M"
        fi
        for x in $(seq 1 ${COR_CONTIG_LAFILTER_RESOLVE_REPEATS})
        do
            tmp="${tmp}${mode}"
        done

        if [[ -z ${FIX_COV} ]]
        then 
            (>&2 echo "If COR_CONTIG_LAFILTER_RESOLVE_REPEATS is set, then the FIX_COV variable has to be set with the apropriate coverage of the patched database ${FIX_DB%.db}.db!")
            exit 1
        fi 
        CONTIG_LAFILTER_OPT="${CONTIG_LAFILTER_OPT} -${tmp} ${FIX_COV}"
    fi    
    if [[ -n ${COR_CONTIG_LAFILTER_STITCH} && ${COR_CONTIG_LAFILTER_STITCH} -gt 0 ]]
    then
        if [[ -n ${COR_CONTIG_LAFILTER_STITCH_AGG} && ${COR_CONTIG_LAFILTER_STITCH_AGG} -gt 0 ]]
        then
            CONTIG_LAFILTER_OPT="${CONTIG_LAFILTER_OPT} -S ${COR_CONTIG_LAFILTER_STITCH}"
        else
            CONTIG_LAFILTER_OPT="${CONTIG_LAFILTER_OPT} -s ${COR_CONTIG_LAFILTER_STITCH}"
        fi
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

function setLAseparateOptions()
{
    CONTIG_LASEPARATE_OPT=""
    if [[ -n ${COR_CONTIG_LASEPARATE_OLEN} && ${COR_CONTIG_LASEPARATE_OLEN} -gt 0 ]]
    then 
        CONTIG_LASEPARATE_OPT="${CONTIG_LASEPARATE_OPT} -o${COR_CONTIG_LASEPARATE_OLEN}"
    fi
    if [[ -n ${COR_CONTIG_LASEPARATE_RLEN} && ${COR_CONTIG_LASEPARATE_RLEN} -gt 0 ]]
    then 
        CONTIG_LASEPARATE_OPT="${CONTIG_LASEPARATE_OPT} -l${COR_CONTIG_LASEPARATE_RLEN}"
    fi 
    if [[ -n ${COR_CONTIG_LASEPARATE_REPEAT} ]]
    then 
        CONTIG_LASEPARATE_OPT="${CONTIG_LASEPARATE_OPT} -r${COR_CONTIG_LASEPARATE_REPEAT}"
    fi 

    # type is passed as argument
    CONTIG_LASEPARATE_OPT="${CONTIG_LASEPARATE_OPT} -T$1"
}

function setForcealignOptions()
{
    CONTIG_FORCEALIGN_OPT=""
    if [[ -n ${COR_CONTIG_FORCEALIGN_PARTIAL} && ${COR_CONTIG_FORCEALIGN_PARTIAL} -ne 0 ]]
    then
        CONTIG_FORCEALIGN_OPT="${CONTIG_FORCEALIGN_OPT} --partial"
    fi
    if [[ -n ${COR_CONTIG_FORCEALIGN_THREADS} && ${COR_CONTIG_FORCEALIGN_THREADS} -gt 0 ]]
    then 
        CONTIG_FORCEALIGN_OPT="${CONTIG_FORCEALIGN_OPT} -t${COR_CONTIG_FORCEALIGN_THREADS}"
    fi 

    if [[ -n ${COR_CONTIG_FORCEALIGN_MAXDIST} && ${COR_CONTIG_FORCEALIGN_MAXDIST} -gt 0 ]]
    then 
        CONTIG_FORCEALIGN_OPT="${CONTIG_FORCEALIGN_OPT} --maxdist${COR_CONTIG_FORCEALIGN_MAXDIST}"
    fi 
    if [[ -n ${COR_CONTIG_FORCEALIGN_BORDER} && ${COR_CONTIG_FORCEALIGN_BORDER} -gt 0 ]]
    then 
        CONTIG_FORCEALIGN_OPT="${CONTIG_FORCEALIGN_OPT} --border${COR_CONTIG_FORCEALIGN_BORDER}"
    fi 
    if [[ -n ${COR_CONTIG_FORCEALIGN_CORRELATION} ]]
    then 
        CONTIG_FORCEALIGN_OPT="${CONTIG_FORCEALIGN_OPT} --correlation${COR_CONTIG_FORCEALIGN_CORRELATION}"
    fi 


    if [[ -z ${COR_CONTIG_FORCEALIGN_RUNID} || ${COR_CONTIG_FORCEALIGN_RUNID} -eq ${COR_CONTIG_DALIGNER_RUNID} ]]
    then
        if [[ -z ${COR_CONTIG_DALIGNER_RUNID} ]]
        then
            COR_CONTIG_FORCEALIGN_RUNID=2
        else 
            COR_CONTIG_FORCEALIGN_RUNID=$((${COR_CONTIG_DALIGNER_RUNID}+1))
        fi
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
    
    ptype="_forcealign"

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
        tmp="${tmp} -t repeats_c${COR_CONTIG_LAREPEAT_COV[$x]}_l${COR_CONTIG_LAREPEAT_LEAVE_COV[$x]}h${COR_CONTIG_LAREPEAT_ENTER_COV[$x]}${ptype}"
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
    if [[ -n ${COR_CONTIG_CTANALYZE_MINCLEN} && ${COR_CONTIG_CTANALYZE_MINCLEN} -gt 0 ]]
    then
        CONTIG_CTANALYZE_OPT="${CONTIG_CTANALYZE_OPT} -x ${COR_CONTIG_CTANALYZE_MINCLEN}"
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
    if [[ -n ${COR_CONTIG_CTANALYZE_REPEATTRACK} && ${COR_CONTIG_CTANALYZE_REPEATTRACK} -gt 0 ]]
    then
        CONTIG_CTANALYZE_OPT="${CONTIG_CTANALYZE_OPT} -r ${COR_CONTIG_CTANALYZE_REPEATTRACK}"
    fi	
    if [[ -n ${COR_CONTIG_CTANALYZE_TRIMTRACK_PATCHEDREADS} && ${COR_CONTIG_CTANALYZE_TRIMTRACK_PATCHEDREADS} -gt 0 ]]
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

myTypes=("1-createCorrectedContigDB, 2-DBdust, 3-Catrack, 4-datander, 5-TANmask, 6-Catrack, 7-daligner, 8-LAfilter, 9-LAseparate, 10-forcealign, 11-LAmerge, 12-LArepeat, 13-TKmerge, 14-TKcombine, 15-LAfilter, 16-LAmerge, 17-CTanalyze")
#type-0 steps: 1-createCorrectedContigDB, 2-DBdust, 3-Catrack, 4-datander, 5-TANmask, 6-Catrack, 7-daligner, 8-LAfilter, 9-LAseparate, 10-forcealign, 11-LAmerge, 12-LArepeat, 13-TKmerge, 14-TKcombine, 15-LAfilter, 16-LAmerge, 17-CTanalyze
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
            if [[ -n ${COR_CONTIG_DALIGNER_NUMACTL} && ${COR_CONTIG_DALIGNER_NUMACTL} -gt 0 ]]
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
                    if [[ -n ${COR_CONTIG_DALIGNER_NUMACTL} && ${COR_CONTIG_DALIGNER_NUMACTL} -gt 0 ]]
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
	### LAfilter - identity overlaps   
    elif [[ ${currentStep} -eq 8 ]]
    then
        ### clean up plans 
        for x in $(ls cont_08_*_*_${CONT_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 

        mkdir -p ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/identity
    	### create LAfilter commands - filter out identity overlaps - has to be done because repcomp and forcealign will loose those
    	for x in $(seq 1 ${contigblocks})
        do  
            echo "${MARVEL_PATH}/bin/LAfilter -p -R 3 ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${CONT_DB%.db} ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/$(getSubDirName ${COR_CONTIG_DALIGNER_RUNID} ${x})/${CONT_DB%.db}.${x}.${CONT_DB%.db}.${x}.las ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/identity/${CONT_DB%.db}.${x}.identity.las"
		done > cont_08_LAfilter_block_${CONT_DB%.db}.${slurmID}.plan
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > cont_08_LAfilter_block_${CONT_DB%.db}.${slurmID}.version
	#### LAseparate 
    elif [[ ${currentStep} -eq 9 ]]
    then
        ### clean up plans 
        for x in $(ls cont_09_*_*_${CONT_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 

        ### find and set forcealign options         
        setDalignerOptions
        setLAseparateOptions 1
        for x in $(seq 1 ${contigblocks}); 
        do 
            sdir=${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/$(getSubDirName ${COR_CONTIG_DALIGNER_RUNID} ${x})
            mkdir -p ${sdir}_ForForceAlign
            mkdir -p ${sdir}_NoForceAlign
            for y in $(seq 1 ${contigblocks}); 
            do 
                if [[ ! -f ${sdir}/${CONT_DB%.db}.${x}.${CONT_DB%.db}.${y}.las ]]
                then
                    (>&2 echo "step ${currentStep} in COR_CONTIG_TYPE ${COR_CONTIG_TYPE}: File missing ${sdir}/${CONT_DB%.db}.${x}.${CONT_DB%.db}.${y}.las!!")
                    exit 1                    
                fi
                echo "${MARVEL_PATH}/bin/LAseparate${CONTIG_LASEPARATE_OPT} ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${CONT_DB%.db} ${sdir}/${CONT_DB%.db}.${x}.${CONT_DB%.db}.${y}.las ${sdir}_ForForceAlign/${CONT_DB%.db}.${x}.${CONT_DB%.db}.${y}.las ${sdir}_NoForceAlign/${CONT_DB%.db}.${x}.${CONT_DB%.db}.${y}.las"                
            done 
		done > cont_09_LAseparate_block_${CONT_DB%.db}.${slurmID}.plan
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > cont_09_LAseparate_block_${CONT_DB%.db}.${slurmID}.version
    #### forcealign 
    elif [[ ${currentStep} -eq 10 ]]
    then
        ### clean up plans 
        for x in $(ls cont_10_*_*_${CONT_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        
        ### find and set forcealign options 
        setForcealignOptions
        cmdLine=1
        for x in $(seq 1 ${contigblocks}); 
        do 
            srcDir=${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/$(getSubDirName ${COR_CONTIG_DALIGNER_RUNID} ${x})_ForForceAlign
            desDir=${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/$(getSubDirName ${COR_CONTIG_FORCEALIGN_RUNID} ${x})

            if [[ ! -d ${desDir} ]]
            then
                mkdir -p ${desDir}
            fi
            start=${x}

            for y in $(seq ${start} ${contigblocks}); 
            do 
                movDir=${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/$(getSubDirName ${COR_CONTIG_FORCEALIGN_RUNID} ${y})
                inFile=${srcDir}/${CONT_DB%.db}.${x}.${CONT_DB%.db}.${y}.las
                
            if [[ -f ${inFile} ]]
                then 
                    if [[ -n ${COR_CONTIG_FORCEALIGN_NUMACTL} && ${COR_CONTIG_FORCEALIGN_NUMACTL} -gt 0 ]]
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
                    echo -n "${NUMACTL}${DACCORD_PATH}/bin/forcealign${CONTIG_FORCEALIGN_OPT} -T/tmp/${CONT_DB%.db}.forcealign.${x}.${y} ${desDir}/${CONT_DB%.db}.forcealign.${x}.${y} ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${CONT_DAZZ_DB%.db} ${inFile}"
                    cmdLine=$((${cmdLine}+1))
                    if [[ $x -eq $y ]]
                    then
                        echo ""
                    else    
                    echo " && mv ${desDir}/${CONT_DB%.db}.forcealign.${x}.${y}_r.las ${movDir}/"
                    fi
                else
                    (>&2 echo "step ${currentStep} in COR_CONTIG_TYPE ${COR_CONTIG_TYPE}: File missing ${inFile}!!")
                    exit 1
                fi
            done 
		done > cont_10_forcealign_block_${CONT_DB%.db}.${slurmID}.plan
		echo "forcealign $(git --git-dir=${DACCORD_SOURCE_PATH}/.git rev-parse --short HEAD)" > cont_10_forcealign_block_${CONT_DB%.db}.${slurmID}.version
    #### LAmerge
    elif [[ ${currentStep} -eq 11 ]]
    then
        ### clean up plans 
        for x in $(ls cont_11_*_*_${CONT_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        
        ### find and set LAmerge options 
        setLAmergeOptions
        setForcealignOptions
        ### create LAmerge commands
        for x in $(seq 1 ${contigblocks}); 
        do  
            echo "${MARVEL_PATH}/bin/LAmerge${CONTIG_LAMERGE_OPT} ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${CONT_DB%.db} ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${CONT_DB%.db}.${x}.forcealign.las ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/$(getSubDirName ${COR_CONTIG_FORCEALIGN_RUNID} ${x}) ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/$(getSubDirName ${COR_CONTIG_DALIGNER_RUNID} ${x})_NoForceAlign ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/identity"
		done > cont_11_LAmerge_block_${CONT_DB%.db}.${slurmID}.plan
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > cont_11_LAmerge_block_${CONT_DB%.db}.${slurmID}.version   
    #### LArepeat
    elif [[ ${currentStep} -eq 12 ]]
    then
        ### clean up plans 
        for x in $(ls cont_12_*_*_${CONT_DB%.db}.${slurmID}.* 2> /dev/null)
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
                echo "${MARVEL_PATH}/bin/LArepeat${CONTIG_LAREPEAT_OPT[$x]} -b ${y} ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${CONT_DB%.db} ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${CONT_DB%.db}.${y}.forcealign.las"
            done
		done > cont_12_LArepeat_block_${CONT_DB%.db}.${slurmID}.plan
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > cont_12_LArepeat_block_${CONT_DB%.db}.${slurmID}.version      
    #### TKmerge 
    elif [[ ${currentStep} -eq 13 ]]
    then
        ### clean up plans 
        for x in $(ls cont_13_*_*_${CONT_DB%.db}.${slurmID}.* 2> /dev/null)
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
		done > cont_13_TKmerge_single_${CONT_DB%.db}.${slurmID}.plan  
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > cont_13_TKmerge_single_${CONT_DB%.db}.${slurmID}.version     
    #### TKcombine 
    elif [[ ${currentStep} -eq 14 ]]
    then
        ### clean up plans 
        for x in $(ls cont_14_*_*_${CONT_DB%.db}.${slurmID}.* 2> /dev/null)
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
		done > cont_14_TKcombine_single_${CONT_DB%.db}.${slurmID}.plan
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > cont_14_TKcombine_single_${CONT_DB%.db}.${slurmID}.version         
	### LAfilter
    elif [[ ${currentStep} -eq 15 ]]
    then
        ### clean up plans 
        for x in $(ls cont_15_*_*_${CONT_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
 
        if [[ -z ${CONTIG_LAFILTER_OPT} ]]
        then 
            setLAfilterOptions
        fi

        ### create LAfilter commands
        contigblocks=$(getNumOfDbBlocks ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${CONT_DB%.db}.db)       
        for x in $(seq 1 ${contigblocks});
        do 
            echo "${MARVEL_PATH}/bin/LAfilter${CONTIG_LAFILTER_OPT} ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${CONT_DB%.db} ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${CONT_DB%.db}.${x}.forcealign.las ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${CONT_DB%.db}.${x}.filt.las"
		done > cont_15_LAfilter_block_${CONT_DB%.db}.${slurmID}.plan 
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > cont_15_LAfilter_block_${CONT_DB%.db}.${slurmID}.version
	### LAmerge
    elif [[ ${currentStep} -eq 16 ]]
    then
        ### clean up plans 
        for x in $(ls cont_16_*_*_${CONT_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        
        ### find and set LAmerge options 
        setLAmergeOptions
        
        echo "${MARVEL_PATH}/bin/LAmerge${CONT_LAMERGE_OPT} -S filt ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${CONT_DB%.db} ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${CONT_DB%.db}.filt.las" > cont_16_LAmerge_single_${CONT_DB%.db}.${slurmID}.plan
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > cont_16_LAmerge_single_${CONT_DB%.db}.${slurmID}.version
	### CTanalyze
	elif [[ ${currentStep} -eq 17 ]]
    then
		### clean up plans 
        for x in $(ls cont_17_*_*_${CONT_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        
        ### find and set CTanalyze options 
        setCTanalyzeOptions
		echo "${MARVEL_PATH}/bin/CTanalyze${CONT_CTANALYZE_OPT} -C ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${CONT_DB%.db} ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${CONT_DB%.db}.filt.las -F ${FIX_FILT_OUTDIR}/${FIX_DB%.DB} ${FIX_FILT_OUTDIR}/${FIX_DB%.DB}.filt.las -D ${FIX_FILT_OUTDIR}/${COR_DIR}/${COR_DB%.db}" > cont_17_CTanalyze_single_${CONT_DB%.db}.${slurmID}.plan
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > cont_17_CTanalyze_single_${CONT_DB%.db}.${slurmID}.version
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