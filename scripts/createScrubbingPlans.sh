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

if [[ ! -n "${FIX_SCRUB_TYPE}" ]]
then 
    (>&2 echo "cannot create read patching scripts if variable FIX_SCRUB_TYPE is not set.")
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

if [[ ${currentStep} -gt 1 && ! -f ${FIX_DB%.db}.db ]]
then 
    (>&2 echo "patched database ${FIX_DB%.db}.db missing")
    exit 1
fi

if [[ ${currentStep} -gt 1 && ${FIX_SCRUB_TYPE} -gt 1 && ! -f ${FIX_DAZZ_DB%.db}.db ]]
then 
    (>&2 echo "patched dazzler database ${RAW_DB%.db}.db missing")
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

rawblocks=$(getNumOfDbBlocks ${RAW_DB%.db}.db)

function setLAfixOptions()
{
    ptype=""
    if [[ ${RAW_PATCH_TYPE} -eq 0 ]]
    then 
        ptype="_dalign"
    elif [[ ${RAW_PATCH_TYPE} -eq 1 ]]
    then 
        ptype="_repcomp"
    elif [[ ${RAW_PATCH_TYPE} -eq 2 ]]
    then 
        ptype="_forcealign"        
    else
        (>&2 echo "Unknown RAW_PATCH_TYPE=${RAW_PATCH_TYPE} !!!")
        exit 1            
    fi 

    if [[ -z ${RAW_FIX_LAFIX_PATH} ]]
    then
        RAW_FIX_LAFIX_PATH=patchedReads${ptype}
    fi 
}

function setdatanderOptions()
{
    ### find and set datander options 
    SCRUB_DATANDER_OPT=""
    if [[ -n ${FIX_SCRUB_DATANDER_THREADS} ]]
    then
        SCRUB_DATANDER_OPT="${SCRUB_DATANDER_OPT} -j ${FIX_SCRUB_DATANDER_THREADS}"
    fi
    if [[ -n ${FIX_SCRUB_DATANDER_MINLEN} ]]
    then
        SCRUB_DATANDER_OPT="${SCRUB_DATANDER_OPT} -l ${FIX_SCRUB_DATANDER_MINLEN}"
    fi
    if [[ -n ${FIX_SCRUB_DATANDER_FOLDER} ]]
    then
        SCRUB_DATANDER_OPT="${SCRUB_DATANDER_OPT} -o ${FIX_SCRUB_DATANDER_FOLDER}"
    else
        FIX_SCRUB_DATANDER_FOLDER="tan"
        SCRUB_DATANDER_OPT="${SCRUB_DATANDER_OPT} -o ${FIX_SCRUB_DATANDER_FOLDER}"
    fi
}

function setTANmaskOptions()
{
    ### find and set TANmask options 
    SCRUB_TANMASK_OPT=""
    if [[ -n ${FIX_SCRUB_TANMASK_VERBOSE} && ${FIX_SCRUB_TANMASK_VERBOSE} -ge 1 ]]
    then
        SCRUB_TANMASK_OPT="${SCRUB_TANMASK_OPT} -v"
    fi
    if [[ -n ${FIX_SCRUB_TANMASK_MINLEN} && ${FIX_SCRUB_TANMASK_MINLEN} -ge 1 ]]
    then
        SCRUB_TANMASK_OPT="${SCRUB_TANMASK_OPT} -l ${FIX_SCRUB_TANMASK_MINLEN}"
    fi
    if [[ -n ${FIX_REPMASK_TANMASK_TRACK} ]]
    then
        SCRUB_TANMASK_OPT="${SCRUB_TANMASK_OPT} -m ${FIX_REPMASK_TANMASK_TRACK}"
    else
        FIX_REPMASK_TANMASK_TRACK="tan"    
        SCRUB_TANMASK_OPT="${SCRUB_TANMASK_OPT} -m ${FIX_REPMASK_TANMASK_TRACK}"
    fi
}

function setCatrackOptions()
{
    ### find and set Catrack options 
    SCRUB_CATRACK_OPT=""
    if [[ -n ${FIX_SCRUB_CATRACK_VERBOSE} && ${FIX_SCRUB_CATRACK_VERBOSE} -ge 1 ]]
    then
        SCRUB_CATRACK_OPT="${SCRUB_CATRACK_OPT} -v"
    fi
    if [[ -n ${FIX_SCRUB_CATRACK_DELETE} && ${FIX_SCRUB_CATRACK_DELETE} -ge 1 ]]
    then
        SCRUB_CATRACK_OPT="${SCRUB_CATRACK_OPT} -d"
    fi
    if [[ -n ${FIX_SCRUB_CATRACK_OVERWRITE} && ${FIX_SCRUB_CATRACK_OVERWRITE} -ge 1 ]]
    then
        SCRUB_CATRACK_OPT="${SCRUB_CATRACK_OPT} -f"
    fi
}

function setDBdustOptions
{
    ### find and set DBdust options 
    SCRUB_DBDUST_OPT=""
    if [[ -n ${FIX_SCRUB_DBDUST_BIAS} && ${FIX_SCRUB_DBDUST_BIAS} -ge 1 ]]
    then
        SCRUB_DBDUST_OPT="${SCRUB_DBDUST_OPT} -b"
    fi
}

function setDalignerOptions()
{
    SCRUB_DALIGNER_OPT=""
    if [[ -n ${FIX_SCRUB_DALIGNER_IDENTITY_OVLS} && ${FIX_SCRUB_DALIGNER_IDENTITY_OVLS} -gt 0 ]]
    then
        SCRUB_DALIGNER_OPT="${SCRUB_DALIGNER_OPT} -I"
    fi
    if [[ -n ${FIX_SCRUB_DALIGNER_KMER} && ${FIX_SCRUB_DALIGNER_KMER} -gt 0 ]]
    then
        SCRUB_DALIGNER_OPT="${SCRUB_DALIGNER_OPT} -k ${FIX_SCRUB_DALIGNER_KMER}"
    fi
    if [[ -n ${FIX_SCRUB_DALIGNER_ERR} ]]
    then
        SCRUB_DALIGNER_OPT="${SCRUB_DALIGNER_OPT} -e ${FIX_SCRUB_DALIGNER_ERR}"
    fi
    if [[ -n ${FIX_SCRUB_DALIGNER_BIAS} && ${FIX_SCRUB_DALIGNER_BIAS} -eq 1 ]]
    then
        SCRUB_DALIGNER_OPT="${SCRUB_DALIGNER_OPT} -b"
    fi
    if [[ -n ${FIX_SCRUB_DALIGNER_VERBOSE} && ${FIX_SCRUB_DALIGNER_VERBOSE} -ne 0 ]]
    then
        SCRUB_DALIGNER_OPT="${SCRUB_DALIGNER_OPT} -v"
    fi
    if [[ -n ${FIX_SCRUB_DALIGNER_TRACESPACE} && ${FIX_SCRUB_DALIGNER_TRACESPACE} -gt 0 ]]
    then
        SCRUB_DALIGNER_OPT="${SCRUB_DALIGNER_OPT} -s ${FIX_SCRUB_DALIGNER_TRACESPACE}"
    fi
    if [[ -n ${FIX_SCRUB_DALIGNER_RUNID} ]]
    then
        SCRUB_DALIGNER_OPT="${SCRUB_DALIGNER_OPT} -r ${FIX_SCRUB_DALIGNER_RUNID}"
    fi
    if [[ -n ${FIX_SCRUB_DALIGNER_T} ]]
    then
        SCRUB_DALIGNER_OPT="${SCRUB_DALIGNER_OPT} -t ${FIX_SCRUB_DALIGNER_T}"
    fi  
    if [[ -n ${FIX_SCRUB_DALIGNER_ASYMMETRIC} && ${FIX_SCRUB_DALIGNER_ASYMMETRIC} -ne 0 ]]
    then
        SCRUB_DALIGNER_OPT="${SCRUB_DALIGNER_OPT} -A"
    fi    
    if [[ -n ${FIX_SCRUB_DALIGNER_MEM} && ${FIX_SCRUB_DALIGNER_MEM} -ne 0 ]]
    then
        SCRUB_DALIGNER_OPT="${SCRUB_DALIGNER_OPT} -M ${FIX_SCRUB_DALIGNER_MEM}"
    fi    
    if [[ -n ${FIX_SCRUB_DALIGNER_MASK} ]]
    then
        for x in ${FIX_SCRUB_DALIGNER_MASK}
        do 
            SCRUB_DALIGNER_OPT="${SCRUB_DALIGNER_OPT} -m ${x}"
        done
    fi
    if [[ -n ${THREADS_daligner} ]]
    then 
        SCRUB_DALIGNER_OPT="${SCRUB_DALIGNER_OPT} -j ${THREADS_daligner}"
    fi
    if [ ! -n ${FIX_SCRUB_DALIGNER_DAL} ]
    then
        FIX_SCRUB_DALIGNER_DAL=8
    fi 
}

function setLAmergeOptions()
{
    SCRUB_LAMERGE_OPT=""
    if [[ -n ${FIX_SCRUB_LAMERGE_NFILES} && ${FIX_SCRUB_LAMERGE_NFILES} -gt 0 ]]
    then
        SCRUB_LAMERGE_OPT="${SCRUB_LAMERGE_OPT} -n ${FIX_SCRUB_LAMERGE_NFILES}"
    fi
}

function setLArepeatOptions()
{
    if [[ ${#FIX_SCRUB_LAREPEAT_LEAVE_COV[*]} -ne ${#FIX_SCRUB_LAREPEAT_ENTER_COV[*]} || ${#FIX_SCRUB_LAREPEAT_ENTER_COV[*]} -ne ${#FIX_SCRUB_LAREPEAT_COV[*]} ]]
    then 
        (>&2 echo "LArepeat number of elements of FIX_SCRUB_LAREPEAT_LEAVE_COV and FIX_SCRUB_LAREPEAT_ENTER_COV and FIX_SCRUB_LAREPEAT_COV differs")
        (>&2 echo "they must be of the same length")
        exit 1
    fi

    numRepeatTracks=${#FIX_SCRUB_LAREPEAT_LEAVE_COV[*]}

    # define array variable - because we may want to create several repeat tracks in one run
    unset SCRUB_LAREPEAT_OPT
    ### find and set LArepeat options     
    
    stype=""
    if [[ "x$1" == "x1" ]]
    then 
        stype="_dalign"
    elif [[ "x$1" == "x2" ]]
    then 
        stype="_repcomp"
    elif [[ "x$1" == "x3" ]]
    then 
        stype="_forcealign"        
    else
        (>&2 echo "Unknown scrubbing type !!!")
        exit 1            
    fi 

    for x in $(seq 0 $((${numRepeatTracks}-1)))
    do 
        tmp=""
        tmp="${tmp} -l ${FIX_SCRUB_LAREPEAT_LEAVE_COV[$x]}"
        tmp="${tmp} -h ${FIX_SCRUB_LAREPEAT_ENTER_COV[$x]}"

        if [[ -n ${FIX_SCRUB_LAREPEAT_OLEN} && ${FIX_SCRUB_LAREPEAT_OLEN} -gt 0 ]]
        then
            tmp="${tmp} -o ${FIX_SCRUB_LAREPEAT_OLEN}"
        fi

        if [[ ${FIX_SCRUB_LAREPEAT_COV[$x]} -ne -1 ]]
        then 
            tmp="${tmp} -c ${FIX_SCRUB_LAREPEAT_COV[$x]}"
            tmp="${tmp} -t repeats_c${FIX_SCRUB_LAREPEAT_COV[$x]}_l${FIX_SCRUB_LAREPEAT_LEAVE_COV[$x]}h${FIX_SCRUB_LAREPEAT_ENTER_COV[$x]}${stype}"
        else
        	if [[ -n ${FIX_SCRUB_LAREPEAT_MAX_COV} && ${FIX_SCRUB_LAREPEAT_MAX_COV} -gt 100 ]]
        	then 
        		tmp="${tmp} -M ${FIX_SCRUB_LAREPEAT_MAX_COV}"
			elif [[ -n ${RAW_COV} && $((${RAW_COV}+20)) -gt 100 ]]
			then
				tmp="${tmp} -M 200"
        	fi         	
            tmp="${tmp} -t repeats_calCov_l${FIX_SCRUB_LAREPEAT_LEAVE_COV[$x]}h${FIX_SCRUB_LAREPEAT_ENTER_COV[$x]}${stype}"
        fi
        SCRUB_LAREPEAT_OPT[$x]=${tmp}
    done 

    FIX_REPMASK_REPEATTRACK=""
    for x in $(seq 1 ${#FIX_REPMASK_BLOCKCMP[*]})
    do
        idx=$(($x-1))
        FIX_REPMASK_REPEATTRACK="${FIX_REPMASK_REPEATTRACK} ${FIX_REPMASK_LAREPEAT_REPEATTRACK}_B${FIX_REPMASK_BLOCKCMP[${idx}]}C${FIX_REPMASK_LAREPEAT_COV[${idx}]}"
    done 

    ## check if repmaskFull_B10C10 exists 
    if [[ -f .${FIX_DB}.${FIX_REPMASK_LAREPEAT_REPEATTRACK}Full_B${FIX_REPMASK_BLOCKCMP[${idx}]}C${FIX_REPMASK_LAREPEAT_COV[${idx}]}.d2 ]]
    then
        FIX_REPMASK_REPEATTRACK="${FIX_REPMASK_REPEATTRACK} ${FIX_REPMASK_LAREPEAT_REPEATTRACK}Full_B${FIX_REPMASK_BLOCKCMP[${idx}]}C${FIX_REPMASK_LAREPEAT_COV[${idx}]}"
    fi    
}

function setTKmergeOptions() 
{
    SCRUB_TKMERGE_OPT=""
    if [[ -n ${FIX_SCRUB_TKMERGE_DELETE} && ${FIX_SCRUB_TKMERGE_DELETE} -ne 0 ]]
    then
        SCRUB_TKMERGE_OPT="${SCRUB_TKMERGE_OPT} -d"
    fi
}

function setTKcombineOptions() 
{
    ignoreDelete=$1
    SCRUB_TKCOMBINE_OPT=""
    if [[ ${ignoreDelete} -eq 0 && -n ${FIX_SCRUB_TKCOMBINE_DELETE} && ${FIX_SCRUB_TKCOMBINE_DELETE} -ne 0 ]]
    then
        SCRUB_TKCOMBINE_OPT="${SCRUB_TKCOMBINE_OPT} -d"
    fi
    if [[ -n ${FIX_SCRUB_TKCOMBINE_VERBOSE} && ${FIX_SCRUB_TKCOMBINE_VERBOSE} -ne 0 ]]
    then
        SCRUB_TKCOMBINE_OPT="${SCRUB_TKCOMBINE_OPT} -v"
    fi
}

function setLAstitchOptions()
{   
    SCRUB_STITCH_OPT=""
    if [[ -n ${FIX_SCRUB_LASTITCH_FUZZ} && ${FIX_SCRUB_LASTITCH_FUZZ} -ne 0 ]]
    then
        SCRUB_STITCH_OPT="${SCRUB_STITCH_OPT} -f ${FIX_SCRUB_LASTITCH_FUZZ}"
    fi
    if [[ -n ${FIX_SCRUB_LASTITCH_ANCHOR} && ${FIX_SCRUB_LASTITCH_ANCHOR} -ne 0 ]]
    then
        SCRUB_STITCH_OPT="${SCRUB_STITCH_OPT} -a ${FIX_SCRUB_LASTITCH_ANCHOR}"
    fi
    if [[ -n ${FIX_SCRUB_LASTITCH_PRELOAD} && ${FIX_SCRUB_LASTITCH_PRELOAD} -ne 0 ]]
    then
        SCRUB_STITCH_OPT="${SCRUB_STITCH_OPT} -L"
    fi
    if [[ -n ${FIX_SCRUB_LASTITCH_PURGE} && ${FIX_SCRUB_LASTITCH_PURGE} -ne 0 ]]
    then
        SCRUB_STITCH_OPT="${SCRUB_STITCH_OPT} -p"
    fi
    if [[ -n ${FIX_SCRUB_LASTITCH_LOWCOMPLEXITY} ]]
    then
        SCRUB_STITCH_OPT="${SCRUB_STITCH_OPT} -t ${FIX_SCRUB_LASTITCH_LOWCOMPLEXITY}"
    fi
    if [[ -n ${FIX_SCRUB_LASTITCH_VERBOSE} ]]
    then
        SCRUB_STITCH_OPT="${SCRUB_STITCH_OPT} -v"
    fi
    if [[ -n ${FIX_SCRUB_LASTITCH_MERGE_DISTANCE} && ${FIX_SCRUB_LASTITCH_MERGE_DISTANCE} -ne 0 ]]
    then
        SCRUB_STITCH_OPT="${SCRUB_STITCH_OPT} -m ${FIX_SCRUB_LASTITCH_MERGE_DISTANCE}"
    fi
    if [[ -n ${FIX_SCRUB_LASTITCH_MINMERGELEN} && ${FIX_SCRUB_LASTITCH_MINMERGELEN} -ne 0 ]]
    then
        SCRUB_STITCH_OPT="${SCRUB_STITCH_OPT} -M ${FIX_SCRUB_LASTITCH_MINMERGELEN}"
    fi    

    # we need the name of the repeat track, especially if the plan starts with step9
	setLArepeatOptions $1

    if [[ -n ${FIX_SCRUB_LASTITCH_REPEATIDX} ]]
    then
        if [[ ${numRepeatTracks} -eq 0 || $((${FIX_SCRUB_LASTITCH_REPEATIDX}+1)) -gt ${#SCRUB_LAREPEAT_OPT[*]} ]]
        then 
            exit 1
        fi

        SCRUB_STITCH_OPT="${SCRUB_STITCH_OPT} -r f$(echo ${SCRUB_LAREPEAT_OPT[${FIX_SCRUB_LASTITCH_REPEATIDX}]} | awk '{print $NF}')_${FIX_REPMASK_LAREPEAT_REPEATTRACK}_${FIX_REPMASK_TANMASK_TRACK}_dust"
    fi
}

function setLAqOptions()
{
    SCRUB_LAQ_OPT=""
    adaptQTRIMCUTOFF=""    

    if [[ -n ${FIX_SCRUB_LAQ_MINSEG} && ${FIX_SCRUB_LAQ_MINSEG} -ne 0 ]]
    then
        SCRUB_LAQ_OPT="${SCRUB_LAQ_OPT} -s ${FIX_SCRUB_LAQ_MINSEG}"
    else 
        FIX_SCRUB_LAQ_MINSEG=25
        SCRUB_LAQ_OPT="${SCRUB_LAQ_OPT} -s ${FIX_SCRUB_LAQ_MINSEG}"
    fi

    if [[ -n ${FIX_SCRUB_LAQ_QTRIMCUTOFF} && ${FIX_SCRUB_LAQ_QTRIMCUTOFF} -ne 0 ]]
    then
        if [[ -n ${RAW_FIX_DALIGNER_TRACESPACE} && ${RAW_FIX_DALIGNER_TRACESPACE} -ne 100 ]]
        then 
            adaptQTRIMCUTOFF=$(echo "${FIX_SCRUB_LAQ_QTRIMCUTOFF}*${RAW_FIX_DALIGNER_TRACESPACE}/100+1" | bc)
            SCRUB_LAQ_OPT="${SCRUB_LAQ_OPT} -d ${adaptQTRIMCUTOFF}"
        else
            adaptQTRIMCUTOFF=${FIX_SCRUB_LAQ_QTRIMCUTOFF}
            SCRUB_LAQ_OPT="${SCRUB_LAQ_OPT} -d ${adaptQTRIMCUTOFF}"            
        fi
    else 
        if [[ -n ${RAW_FIX_DALIGNER_TRACESPACE} && ${RAW_FIX_DALIGNER_TRACESPACE} -ne 100 ]]
        then 
            FIX_SCRUB_LAQ_QTRIMCUTOFF=25
            adaptQTRIMCUTOFF=$(echo "${FIX_SCRUB_LAQ_QTRIMCUTOFF}*${RAW_FIX_DALIGNER_TRACESPACE}/100+1" | bc)
            SCRUB_LAQ_OPT="${SCRUB_LAQ_OPT} -d ${adaptQTRIMCUTOFF}"
        else
            adaptQTRIMCUTOFF=25
            FIX_SCRUB_LAQ_QTRIMCUTOFF=25
            SCRUB_LAQ_OPT="${SCRUB_LAQ_OPT} -d ${adaptQTRIMCUTOFF}"            
        fi
    fi
}

function setLAgapOptions()
{
    SCRUB_LAGAP_OPT=""
    if [[ -n ${FIX_SCRUB_LAGAP_STITCH} && ${FIX_SCRUB_LAGAP_STITCH} -gt 0 ]]
    then
        SCRUB_LAGAP_OPT="${SCRUB_LAGAP_OPT} -s ${FIX_SCRUB_LAGAP_STITCH}"
    fi
    if [[ -n ${FIX_SCRUB_LAGAP_PRELOAD} && ${FIX_SCRUB_LAGAP_PRELOAD} -ne 0 ]]
    then
        SCRUB_LAGAP_OPT="${SCRUB_LAGAP_OPT} -L"
    fi
    if [[ -n ${FIX_SCRUB_LAGAP_PURGE} && ${FIX_SCRUB_LAGAP_PURGE} -ne 0 ]]
    then
        SCRUB_LAGAP_OPT="${SCRUB_LAGAP_OPT} -p"
    fi
    
    if [[ -n ${FIX_SCRUB_LAGAP_TRIM} && ${FIX_SCRUB_LAGAP_TRIM} -ne 0 ]]
    then
        if [[ -z ${SCRUB_LAQ_OPT} ]]
        then
            setLAqOptions
        fi
        ### should bet set within the plans
        ###SCRUB_LAGAP_OPT="${SCRUB_LAGAP_OPT} -t trim0_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}"
    fi

    if [[ -n ${FIX_SCRUB_LAGAP_DISCARD_CHIMERS} ]]
    then
            setLArepeatOptions $1            

            if [[ ${numRepeatTracks} -eq 0 || $((${FIX_SCRUB_LAGAP_DISCARD_CHIMERS}+1)) -gt ${#SCRUB_LAREPEAT_OPT[*]} ]]
            then 
                exit 1
            fi
            tmp=$(echo ${SCRUB_LAREPEAT_OPT[${FIX_SCRUB_LASTITCH_REPEATIDX}]} | awk '{print $NF}')_${FIX_REPMASK_LAREPEAT_REPEATTRACK}

            SCRUB_LAGAP_OPT="${SCRUB_LAGAP_OPT} -R f${tmp}_${FIX_REPMASK_TANMASK_TRACK}_dust"
    fi
    
}

function setRepcompOptions()
{
    SCRUB_REPCOMP_OPT=""
    if [[ -n ${FIX_SCRUB_REPCOMP_TRACESPACE} && ${FIX_SCRUB_REPCOMP_TRACESPACE} -gt 0 ]]
    then 
        SCRUB_REPCOMP_OPT="${SCRUB_REPCOMP_OPT} --tspace${FIX_SCRUB_REPCOMP_TRACESPACE}"
    fi
    if [[ -n ${FIX_SCRUB_REPCOMP_INBLOCKSIZE} ]]
    then 
        SCRUB_REPCOMP_OPT="${SCRUB_REPCOMP_OPT} -i${FIX_SCRUB_REPCOMP_INBLOCKSIZE}"
    fi
    if [[ -n ${FIX_SCRUB_REPCOMP_KMER} && ${FIX_SCRUB_REPCOMP_KMER} -gt 0 ]]
    then 
        SCRUB_REPCOMP_OPT="${SCRUB_REPCOMP_OPT} -k${FIX_SCRUB_REPCOMP_KMER}"
    fi
    if [[ -n ${FIX_SCRUB_REPCOMP_MEM} ]]
    then 
        SCRUB_REPCOMP_OPT="${SCRUB_REPCOMP_OPT} -M${FIX_SCRUB_REPCOMP_MEM}"
    fi
    if [[ -n ${FIX_SCRUB_REPCOMP_THREADS} && ${FIX_SCRUB_REPCOMP_THREADS} -gt 0 ]]
    then 
        SCRUB_REPCOMP_OPT="${SCRUB_REPCOMP_OPT} -t${FIX_SCRUB_REPCOMP_THREADS}"
    fi 

    if [[ -n ${FIX_SCRUB_REPCOMP_CORRELATION} ]]
    then 
        SCRUB_REPCOMP_OPT="${SCRUB_REPCOMP_OPT} -e${FIX_SCRUB_REPCOMP_CORRELATION}"
    fi 
    if [[ -n ${FIX_SCRUB_REPCOMP_MASK} ]]
    then
        for x in ${FIX_SCRUB_REPCOMP_MASK}
        do 
            SCRUB_REPCOMP_OPT="${SCRUB_REPCOMP_OPT} -m${x}"
        done
    fi
    if [[ -n ${FIX_SCRUB_REPCOMP_OLEN} && ${FIX_SCRUB_REPCOMP_OLEN} -gt 0 ]]
    then 
        SCRUB_REPCOMP_OPT="${SCRUB_REPCOMP_OPT} -l${FIX_SCRUB_REPCOMP_OLEN}"
    fi

    if [[ -z ${FIX_SCRUB_REPCOMP_RUNID} || ${FIX_SCRUB_REPCOMP_RUNID} -eq ${FIX_SCRUB_DALIGNER_RUNID} ]]
    then
        FIX_SCRUB_REPCOMP_RUNID=$((${FIX_SCRUB_REPCOMP_RUNID}+1))
    fi
}

function setForcealignOptions()
{
    SCRUB_FORCEALIGN_OPT=""
    if [[ -n ${FIX_SCRUB_FORCEALIGN_PARTIAL} && ${FIX_SCRUB_FORCEALIGN_PARTIAL} -ne 0 ]]
    then
        SCRUB_FORCEALIGN_OPT="${SCRUB_FORCEALIGN_OPT} --partial"
    fi
    if [[ -n ${FIX_SCRUB_FORCEALIGN_THREADS} && ${FIX_SCRUB_FORCEALIGN_THREADS} -gt 0 ]]
    then 
        SCRUB_FORCEALIGN_OPT="${SCRUB_FORCEALIGN_OPT} -t${FIX_SCRUB_FORCEALIGN_THREADS}"
    fi 

    if [[ -n ${FIX_SCRUB_FORCEALIGN_MAXDIST} && ${FIX_SCRUB_FORCEALIGN_MAXDIST} -gt 0 ]]
    then 
        SCRUB_FORCEALIGN_OPT="${SCRUB_FORCEALIGN_OPT} --maxdist${FIX_SCRUB_FORCEALIGN_MAXDIST}"
    fi 
    if [[ -n ${FIX_SCRUB_FORCEALIGN_BORDER} && ${FIX_SCRUB_FORCEALIGN_BORDER} -gt 0 ]]
    then 
        SCRUB_FORCEALIGN_OPT="${SCRUB_FORCEALIGN_OPT} --border${FIX_SCRUB_FORCEALIGN_BORDER}"
    fi 
    if [[ -n ${FIX_SCRUB_FORCEALIGN_CORRELATION} ]]
    then 
        SCRUB_FORCEALIGN_OPT="${SCRUB_FORCEALIGN_OPT} --correlation${FIX_SCRUB_FORCEALIGN_CORRELATION}"
    fi 


    if [[ -z ${FIX_SCRUB_FORCEALIGN_RUNID} || ${FIX_SCRUB_FORCEALIGN_RUNID} -eq ${FIX_SCRUB_REPCOMP_RUNID} ]]
    then
        if [[ -z ${FIX_SCRUB_REPCOMP_RUNID} ]]
        then
            FIX_SCRUB_FORCEALIGN_RUNID=$((${FIX_SCRUB_DALIGNER_RUNID}+2))
        else 
            FIX_SCRUB_FORCEALIGN_RUNID=$((${FIX_SCRUB_REPCOMP_RUNID}+1))
        fi
    fi
}

function setLAseparateOptions()
{
	SCRUB_LASEPARATE_OPT=""
    if [[ -n ${FIX_SCRUB_LASEPARATE_OLEN} && ${FIX_SCRUB_LASEPARATE_OLEN} -gt 0 ]]
    then 
        SCRUB_LASEPARATE_OPT="${SCRUB_LASEPARATE_OPT} -o${FIX_SCRUB_LASEPARATE_OLEN}"
    fi
    if [[ -n ${FIX_SCRUB_LASEPARATE_RLEN} && ${FIX_SCRUB_LASEPARATE_RLEN} -gt 0 ]]
    then 
        SCRUB_LASEPARATE_OPT="${SCRUB_LASEPARATE_OPT} -l${FIX_SCRUB_LASEPARATE_RLEN}"
    fi 
    if [[ -n ${FIX_SCRUB_LASEPARATE_USEREPEATIDX} && ${FIX_SCRUB_LASEPARATE_USEREPEATIDX} -ge 0 && ${FIX_SCRUB_LASEPARATE_USEREPEATIDX} -lt ${#FIX_SCRUB_LAREPEAT_LEAVE_COV[*]} ]]
    then 
    	stype=""
    	if [[ "x$1" == "x0" ]]
    	then 
    		setLArepeatOptions 1
    	elif [[ "x$1" == "x1" ]]
    	then 
    		setLArepeatOptions 2
    	fi
    		
		FIX_SCRUB_LASEPARATE_REPEAT="f$(echo ${SCRUB_LAREPEAT_OPT[${FIX_SCRUB_LASEPARATE_USEREPEATIDX}]} | awk '{print $NF}')_${FIX_REPMASK_LAREPEAT_REPEATTRACK}_${FIX_REPMASK_TANMASK_TRACK}_dust"
    	SCRUB_LASEPARATE_OPT="${SCRUB_LASEPARATE_OPT} -r${FIX_SCRUB_LASEPARATE_REPEAT}"
    fi
    
    # type is passed as argument
    SCRUB_LASEPARATE_OPT="${SCRUB_LASEPARATE_OPT} -T$1"
}

## ensure some paths
if [[ -z "${MARVEL_SOURCE_PATH}" || ! -d "${MARVEL_SOURCE_PATH}" ]]
then 
    (>&2 echo "ERROR - You have to set MARVEL_SOURCE_PATH. Used to report git version.")
    exit 1
fi

if [[ -z "${DAZZLER_SOURCE_PATH}" ]]
then 
    (>&2 echo "ERROR - You have to set DAZZLER_SOURCE_PATH. Used to report git version.")
    exit 1
fi

if [[ -z "${REPCOMP_SOURCE_PATH}" ]]
then 
    (>&2 echo "ERROR - You have to set REPCOMP_SOURCE_PATH. Used to report git version.")
    exit 1
fi

if [[ -z "${DACCORD_SOURCE_PATH}" ]]
then 
    (>&2 echo "ERROR - You have to set DACCORD_SOURCE_PATH. Used to report git version.")
    exit 1
fi

fixblocks=$(getNumOfDbBlocks ${FIX_DB%.db}.db)

mytypes=("1-daligner, 2-LAmerge, 3-LArepeat, 4-TKmerge, 5-TKcombine, 6-TKhomogenize, 7-TKcombine, 8-LAstitch, 9-LAq, 10-TKmerge, 11-LAgap, 12-LAq, 13-TKmerge" "1-daligner, 2-LAmerge, 3-LArepeat, 4-TKmerge, 5-TKcombine, 6-TKhomogenize, 7-TKcombine, 8-LAstitch, 9-LAq, 10-TKmerge, 11-LAgap, 12-LAq, 13-TKmerge, 14-LAseparate, 15-repcomp, 16-LAmerge, 17-LArepeat, 18-TKmerge, 19-TKcombine, 20-TKhomogenize, 21-TKcombine, 22-LAstitch, 23-LAq, 24-TKmerge, 25-LAgap, 26-LAq, 27-TKmerge, 28-LAseparate, 29-forcealign, 30-LAmerge, 31-LArepeat, 32-TKmerge, 33-TKcombine, 34-TKhomogenize, 35-TKcombine, 36-LAstitch, 37-LAq, 38-TKmerge, 39-LAgap, 40-LAq, 41-TKmerge")

#type-0 steps [1-13]: 1-daligner, 2-LAmerge, 3-LArepeat, 4-TKmerge, 5-TKcombine, 6-TKhomogenize, 7-TKcombine, 8-LAstitch, 9-LAq, 10-TKmerge, 11-LAgap, 12-LAq, 13-TKmerge          ## old pipeline
if [[ ${FIX_SCRUB_TYPE} -eq 0 ]]
then
    #### daligner 
    if [[ ${currentStep} -eq 1 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_01_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        ### find and set daligner options 
        setDalignerOptions
        cmdLine=1
        ### create daligner commands
        for x in $(seq 1 ${fixblocks})
        do 
            if [[ -n ${FIX_SCRUB_DALIGNER_NUMACTL} && ${FIX_SCRUB_DALIGNER_NUMACTL} -gt 0 ]] && [[ "x${SLURM_NUMACTL}" == "x" || ${SLURM_NUMACTL} -eq 0 ]]
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
            cmd="${NUMACTL}${MARVEL_PATH}/bin/daligner${SCRUB_DALIGNER_OPT} ${FIX_DB%.db}.${x}"
            cmdLine=$((${cmdLine}+1))
            count=0
            for y in $(seq ${x} ${fixblocks})
            do  
                if [[ $count -lt ${FIX_SCRUB_DALIGNER_DAL} ]]
                then
                    cmd="${cmd} ${FIX_DB%.db}.${y}"
                    count=$((${count}+1))
                else    
                    echo "${cmd}"
                    if [[ -n ${FIX_SCRUB_DALIGNER_NUMACTL} && ${FIX_SCRUB_DALIGNER_NUMACTL} -gt 0 ]] && [[ "x${SLURM_NUMACTL}" == "x" || ${SLURM_NUMACTL} -eq 0 ]]
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
                    cmd="${NUMACTL}${MARVEL_PATH}/bin/daligner${SCRUB_DALIGNER_OPT} ${FIX_DB%.db}.${x} ${FIX_DB%.db}.${y}"
                    cmdLine=$((${cmdLine}+1))
                    count=1
                fi
            done
            echo "${cmd}"
    	done > scrub_01_daligner_block_${FIX_DB%.db}.${slurmID}.plan
    	echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_01_daligner_block_${FIX_DB%.db}.${slurmID}.version
    #### LAmerge
    elif [[ ${currentStep} -eq 2 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_02_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        
        ### find and set LAmerge options 
        setLAmergeOptions
        ### create LAmerge commands        
        for x in $(seq 1 ${fixblocks})
        do  
            echo "${MARVEL_PATH}/bin/LAmerge${SCRUB_LAMERGE_OPT} ${FIX_DB%.db} ${FIX_DB%.db}.${x}.dalign.las $(getSubDirName ${FIX_SCRUB_DALIGNER_RUNID} ${x})"
    	done > scrub_02_LAmerge_block_${FIX_DB%.db}.${slurmID}.plan
    	echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_02_LAmerge_block_${FIX_DB%.db}.${slurmID}.version   
    #### LArepeat
    elif [[ ${currentStep} -eq 3 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_03_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 

        setLArepeatOptions 1
        if [[ ${numRepeatTracks} -eq 0 ]]
        then 
            exit 1
        fi    
    
        for x in $(seq 0 $((${numRepeatTracks}-1)))
        do 
            ### create LArepeat commands
            for y in $(seq 1 ${fixblocks})
            do 
                echo "${MARVEL_PATH}/bin/LArepeat${SCRUB_LAREPEAT_OPT[$x]} -b ${y} ${FIX_DB%.db} ${FIX_DB%.db}.${y}.dalign.las"
            done
    	done > scrub_03_LArepeat_block_${FIX_DB%.db}.${slurmID}.plan
    	echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_03_LArepeat_block_${FIX_DB%.db}.${slurmID}.version      
    #### TKmerge 
    elif [[ ${currentStep} -eq 4 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_04_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done         

        setLArepeatOptions 1
        ### find and set TKmerge options 
        if [[ -z ${SCRUB_TKMERGE_OPT} ]]
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
            echo "${MARVEL_PATH}/bin/TKmerge${SCRUB_TKMERGE_OPT} ${FIX_DB%.db} $(echo ${SCRUB_LAREPEAT_OPT[${x}]} | awk '{print $NF}')" 
    	 done > scrub_04_TKmerge_single_${FIX_DB%.db}.${slurmID}.plan
    	 echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_04_TKmerge_single_${FIX_DB%.db}.${slurmID}.version       
    #### TKcombine 
    elif [[ ${currentStep} -eq 5 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_05_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done             
         # we need the name of the repeat track, especially if the plan starts with step11
        setLArepeatOptions 1
        
        ### find and set TKcombine options
        setTKcombineOptions 1

        if [[ ${numRepeatTracks} -eq 0 ]]
        then 
            exit 1
        fi 
        ### create TKmerge command
        if [[ -n ${FIX_REPMASK_REPEATTRACK} ]]
        then 
            for x in $(seq 0 $((${numRepeatTracks}-1)))
            do 
                tmp=$(echo ${SCRUB_LAREPEAT_OPT[${x}]} | awk '{print $NF}')
                echo "${MARVEL_PATH}/bin/TKcombine${SCRUB_TKCOMBINE_OPT} ${FIX_DB%.db} ${tmp}_${FIX_REPMASK_LAREPEAT_REPEATTRACK} ${tmp} ${FIX_REPMASK_REPEATTRACK}" 
        	done > scrub_05_TKcombine_single_${FIX_DB%.db}.${slurmID}.plan         
        else
            for x in $(seq 0 $((${numRepeatTracks}-1)))
            do 
                tmp=$(echo ${SCRUB_LAREPEAT_OPT[${x}]} | awk '{print $NF}')
                echo "ln -s -f .${FIX_DB%.db}.${tmp}.d2 .${FIX_DB%.db}.${tmp}_${FIX_REPMASK_LAREPEAT_REPEATTRACK}.d2"
                echo "ln -s -f .${FIX_DB%.db}.${tmp}.a2 .${FIX_DB%.db}.${tmp}_${FIX_REPMASK_LAREPEAT_REPEATTRACK}.a2"

        	done > scrub_05_TKcombine_single_${FIX_DB%.db}.${slurmID}.plan
        fi 
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_05_TKcombine_single_${FIX_DB%.db}.${slurmID}.version
    #### TKhomogenize
    elif [[ ${currentStep} -eq 6 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_06_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done             
         # we need the name of the repeat track, especially if the plan starts with step6
        setLArepeatOptions 1
        ### create TKhomogenize commands
        for x in $(seq 0 $((${numRepeatTracks}-1)))
        do 
            tmp=$(echo ${SCRUB_LAREPEAT_OPT[${x}]} | awk '{print $NF}')_${FIX_REPMASK_LAREPEAT_REPEATTRACK}
            for y in $(seq 1 ${fixblocks})
            do 
                echo "${MARVEL_PATH}/bin/TKhomogenize -i  ${tmp} -I h${tmp} -b ${y} ${FIX_DB%.db} ${FIX_DB%.db}.${y}.dalign.las"
            done 
    	done > scrub_06_TKhomogenize_block_${FIX_DB%.db}.${slurmID}.plan
    	echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_06_TKhomogenize_block_${FIX_DB%.db}.${slurmID}.version         
    #### TKcombine
    elif [[ ${currentStep} -eq 7 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_07_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done             
         # we need the name of the repeat track, especially if the plan starts with step7
        setLArepeatOptions 1 
        ### find and set TKcombine options
        setTKcombineOptions 0
        ### create TKcombine commands
        for x in $(seq 0 $((${numRepeatTracks}-1)))
        do 
            tmp=$(echo ${SCRUB_LAREPEAT_OPT[${x}]} | awk '{print $NF}')_${FIX_REPMASK_LAREPEAT_REPEATTRACK}
            echo "${MARVEL_PATH}/bin/TKcombine${SCRUB_TKCOMBINE_OPT} ${FIX_DB%.db} h${tmp} \#.h${tmp}"
    	done > scrub_07_TKcombine_single_${FIX_DB%.db}.${slurmID}.plan         
        ### find and set TKcombine options, but ignore the -d flag if it was set
        setTKcombineOptions 1

        echo "${MARVEL_PATH}/bin/TKcombine${SCRUB_TKCOMBINE_OPT} ${FIX_DB%.db} ${FIX_REPMASK_TANMASK_TRACK}_dust ${FIX_REPMASK_TANMASK_TRACK} dust" >> scrub_07_TKcombine_single_${FIX_DB%.db}.${slurmID}.plan
        for x in $(seq 0 $((${numRepeatTracks}-1)))
        do
            tmp=$(echo ${SCRUB_LAREPEAT_OPT[${x}]} | awk '{print $NF}')_${FIX_REPMASK_LAREPEAT_REPEATTRACK}
            echo "${MARVEL_PATH}/bin/TKcombine${SCRUB_TKCOMBINE_OPT} ${FIX_DB%.db} f${tmp} h${tmp} ${tmp}" 
            echo "${MARVEL_PATH}/bin/TKcombine${SCRUB_TKCOMBINE_OPT} ${FIX_DB%.db} f${tmp}_${FIX_REPMASK_TANMASK_TRACK} f${tmp} ${FIX_REPMASK_TANMASK_TRACK}" 
            echo "${MARVEL_PATH}/bin/TKcombine${SCRUB_TKCOMBINE_OPT} ${FIX_DB%.db} f${tmp}_${FIX_REPMASK_TANMASK_TRACK}_dust f${tmp}_${FIX_REPMASK_TANMASK_TRACK} dust"            
		done >> scrub_07_TKcombine_single_${FIX_DB%.db}.${slurmID}.plan
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_07_TKcombine_single_${FIX_DB%.db}.${slurmID}.version      
    #### LAstitch
    elif [[ ${currentStep} -eq 8 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_08_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done             
        ### find and set LAstitch options 
        setLAstitchOptions 1
        ### create LAstitch commands
        for x in $(seq 1 ${fixblocks})
        do 
            echo "${MARVEL_PATH}/bin/LAstitch${SCRUB_STITCH_OPT} ${FIX_DB%.db} ${FIX_DB%.db}.${x}.dalign.las ${FIX_DB%.db}.${x}.dalignStitch.las"
		done > scrub_08_LAstitch_block_${FIX_DB%.db}.${slurmID}.plan
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_08_LAstitch_block_${FIX_DB%.db}.${slurmID}.version                 
    #### LAq
    elif [[ ${currentStep} -eq 9 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_09_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 

        ### find and set LAq options 
        setLAqOptions
        ### create LAq commands
        for x in $(seq 1 ${fixblocks})
        do 
            echo "${MARVEL_PATH}/bin/LAq${SCRUB_LAQ_OPT} -T trim0_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}_dalign -Q q0_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}_dalign -b ${x} ${FIX_DB%.db} ${FIX_DB%.db}.${x}.dalignStitch.las"
    	done > scrub_09_LAq_block_${FIX_DB%.db}.${slurmID}.plan
    	echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_09_LAq_block_${FIX_DB%.db}.${slurmID}.version               
    #### TKmerge    
    elif [[ ${currentStep} -eq 10 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_10_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done         
        # we need the name of the repeat track, especially if the plan starts with step10
        if [[ -z ${SCRUB_LAQ_OPT} ]]
        then 
            setLAqOptions
        fi
        ### find and set TKmerge options 
        if [[ -z ${SCRUB_TKMERGE_OPT} ]]
        then 
            setTKmergeOptions
        fi
        
        ### create TKmerge command
        echo "${MARVEL_PATH}/bin/TKmerge${SCRUB_TKMERGE_OPT} ${FIX_DB%.db} trim0_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}_dalign" > scrub_10_TKmerge_block_${FIX_DB%.db}.${slurmID}.plan 
        echo "${MARVEL_PATH}/bin/TKmerge${SCRUB_TKMERGE_OPT} ${FIX_DB%.db} q0_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}_dalign" >> scrub_10_TKmerge_block_${FIX_DB%.db}.${slurmID}.plan
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_10_TKmerge_block_${FIX_DB%.db}.${slurmID}.version               
    #### LAgap
    elif [[ ${currentStep} -eq 11 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_11_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done         
        setLAgapOptions 1
        ### create LAq commands
        for x in $(seq 1 ${fixblocks})
        do 
            breakChim=""
            if [[ -n ${FIX_SCRUB_LAGAP_DISCARD_CHIMERS} ]]
            then 
                breakChim=" -C ${FIX_DB%.db}.${x}.dalignGap.chimers.txt"
            fi
            echo "${MARVEL_PATH}/bin/LAgap${SCRUB_LAGAP_OPT} -t trim0_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}_dalign${breakChim} ${FIX_DB%.db} ${FIX_DB%.db}.${x}.dalignStitch.las ${FIX_DB%.db}.${x}.dalignGap.las"
    	done > scrub_11_LAgap_block_${FIX_DB%.db}.${slurmID}.plan      
    	echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_11_LAgap_block_${FIX_DB%.db}.${slurmID}.version         
    #### LAq
    elif [[ ${currentStep} -eq 12 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_12_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 

        ### find and set LAq options 
        setLAqOptions
        ### create LAq commands
        for x in $(seq 1 ${fixblocks})
        do 
            echo "${MARVEL_PATH}/bin/LAq${SCRUB_LAQ_OPT} -u -T trim1_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}_dalign -t trim0_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}_dalign -q q0_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}_dalign -b ${x} ${FIX_DB%.db} ${FIX_DB%.db}.${x}.dalignGap.las"
    done > scrub_12_LAq_block_${FIX_DB%.db}.${slurmID}.plan 
    echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_12_LAq_block_${FIX_DB%.db}.${slurmID}.version              
    #### TKmerge    
    elif [[ ${currentStep} -eq 13 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_13_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done         
        # we need the name of the repeat track, especially if the plan starts with step10
        if [[ -z ${SCRUB_LAQ_OPT} ]]
        then 
            setLAqOptions
        fi
        ### find and set TKmerge options 
        if [[ -z ${SCRUB_TKMERGE_OPT} ]]
        then 
            setTKmergeOptions
        fi
        
        ### create TKmerge command
        echo "${MARVEL_PATH}/bin/TKmerge${SCRUB_TKMERGE_OPT} ${FIX_DB%.db} trim1_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}_dalign" > scrub_13_TKmerge_single_${FIX_DB%.db}.${slurmID}.plan 
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_13_TKmerge_single_${FIX_DB%.db}.${slurmID}.version           
    else
        (>&2 echo "step ${currentStep} in FIX_SCRUB_TYPE ${FIX_SCRUB_TYPE} not supported")
        (>&2 echo "valid steps are: ${myTypes[${FIX_SCRUB_TYPE}]}")
        exit 1            
    fi
#type-1 steps [ 1-13 : 1-daligner, 2-LAmerge, 3-LArepeat, 4-TKmerge, 5-TKcombine, 6-TKhomogenize, 7-TKcombine, 8-LAstitch, 9-LAq, 10-TKmerge, 11-LAgap, 12-LAq, 13-TKmerge          ## experimental pipeline
#              14-27 : 14-LAseparate, 15-repcomp, 16-LAmerge, 17-LArepeat, 18-TKmerge, 19-TKcombine, 20-TKhomogenize, 21-TKcombine, 22-LAstitch, 23-LAq, 24-TKmerge, 25-LAgap, 26-LAq, 27-TKmerge
#              28-41]: 28-LAseparate, 29-forcealign, 30-LAmerge, 31-LArepeat, 32-TKmerge, 33-TKcombine, 34-TKhomogenize, 35-TKcombine, 36-LAstitch, 37-LAq, 38-TKmerge, 39-LAgap, 40-LAq, 41-TKmerge
elif [[ ${FIX_SCRUB_TYPE} -eq 1 ]]
then
	#### daligner 
    if [[ ${currentStep} -eq 1 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_01_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        ### find and set daligner options 
        setDalignerOptions
        cmdLine=1
        ### create daligner commands
        for x in $(seq 1 ${fixblocks})
        do 
            if [[ -n ${FIX_SCRUB_DALIGNER_NUMACTL} && ${FIX_SCRUB_DALIGNER_NUMACTL} -gt 0 ]] && [[ "x${SLURM_NUMACTL}" == "x" || ${SLURM_NUMACTL} -eq 0 ]]
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
            cmd="${NUMACTL}${MARVEL_PATH}/bin/daligner${SCRUB_DALIGNER_OPT} ${FIX_DB%.db}.${x}"
            cmdLine=$((${cmdLine}+1))
            count=0
            for y in $(seq ${x} ${fixblocks})
            do  
                if [[ $count -lt ${FIX_SCRUB_DALIGNER_DAL} ]]
                then
                    cmd="${cmd} ${FIX_DB%.db}.${y}"
                    count=$((${count}+1))
                else    
                    echo "${cmd}"
                    if [[ -n ${FIX_SCRUB_DALIGNER_NUMACTL} && ${FIX_SCRUB_DALIGNER_NUMACTL} -gt 0 ]] && [[ "x${SLURM_NUMACTL}" == "x" || ${SLURM_NUMACTL} -eq 0 ]]
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
                    cmd="${NUMACTL}${MARVEL_PATH}/bin/daligner${SCRUB_DALIGNER_OPT} ${FIX_DB%.db}.${x} ${FIX_DB%.db}.${y}"
                    cmdLine=$((${cmdLine}+1))
                    count=1
                fi
            done
            echo "${cmd}"
    	done > scrub_01_daligner_block_${FIX_DB%.db}.${slurmID}.plan
    	echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_01_daligner_block_${FIX_DB%.db}.${slurmID}.version
    #### LAmerge
    elif [[ ${currentStep} -eq 2 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_02_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        
        ### find and set LAmerge options 
        setLAmergeOptions
        ### create LAmerge commands
        for x in $(seq 1 ${fixblocks})
        do  
            echo "${MARVEL_PATH}/bin/LAmerge${SCRUB_LAMERGE_OPT} ${FIX_DB%.db} ${FIX_DB%.db}.${x}.dalign.las $(getSubDirName ${FIX_SCRUB_DALIGNER_RUNID} ${x})"
    	done > scrub_02_LAmerge_block_${FIX_DB%.db}.${slurmID}.plan
    	echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_02_LAmerge_block_${FIX_DB%.db}.${slurmID}.version   
    #### LArepeat
    elif [[ ${currentStep} -eq 3 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_03_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 

        setLArepeatOptions 1
        if [[ ${numRepeatTracks} -eq 0 ]]
        then 
            exit 1
        fi    
    
        for x in $(seq 0 $((${numRepeatTracks}-1)))
        do 
            ### create LArepeat commands
            for y in $(seq 1 ${fixblocks})
            do 
                echo "${MARVEL_PATH}/bin/LArepeat${SCRUB_LAREPEAT_OPT[$x]} -b ${y} ${FIX_DB%.db} ${FIX_DB%.db}.${y}.dalign.las"
            done
    	done > scrub_03_LArepeat_block_${FIX_DB%.db}.${slurmID}.plan
    	echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_03_LArepeat_block_${FIX_DB%.db}.${slurmID}.version      
    #### TKmerge 
    elif [[ ${currentStep} -eq 4 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_04_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done         
        # we need the name of the repeat track, especially if the plan starts with step10
        setLArepeatOptions 1
        ### find and set TKmerge options 
        if [[ -z ${SCRUB_TKMERGE_OPT} ]]
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
            echo "${MARVEL_PATH}/bin/TKmerge${SCRUB_TKMERGE_OPT} ${FIX_DB%.db} $(echo ${SCRUB_LAREPEAT_OPT[${x}]} | awk '{print $NF}')" 
    	 done > scrub_04_TKmerge_single_${FIX_DB%.db}.${slurmID}.plan
    	 echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_04_TKmerge_single_${FIX_DB%.db}.${slurmID}.version       
    #### TKcombine 
    elif [[ ${currentStep} -eq 5 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_05_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done             
         # we need the name of the repeat track, especially if the plan starts with step11
        setLArepeatOptions 1
        ### find and set TKcombine options
        setTKcombineOptions 1

        if [[ ${numRepeatTracks} -eq 0 ]]
        then 
            exit 1
        fi 
        ### create TKmerge command
        if [[ -n ${FIX_REPMASK_REPEATTRACK} ]]
        then 
            for x in $(seq 0 $((${numRepeatTracks}-1)))
            do 
                tmp=$(echo ${SCRUB_LAREPEAT_OPT[${x}]} | awk '{print $NF}')
                echo "${MARVEL_PATH}/bin/TKcombine${SCRUB_TKCOMBINE_OPT} ${FIX_DB%.db} ${tmp}_${FIX_REPMASK_LAREPEAT_REPEATTRACK} ${tmp} ${FIX_REPMASK_REPEATTRACK}" 
        	done > scrub_05_TKcombine_single_${FIX_DB%.db}.${slurmID}.plan         
        else
            for x in $(seq 0 $((${numRepeatTracks}-1)))
            do 
                tmp=$(echo ${SCRUB_LAREPEAT_OPT[${x}]} | awk '{print $NF}')
                echo "ln -s -f .${FIX_DB%.db}.${tmp}.d2 .${FIX_DB%.db}.${tmp}_${FIX_REPMASK_LAREPEAT_REPEATTRACK}.d2"
                echo "ln -s -f .${FIX_DB%.db}.${tmp}.a2 .${FIX_DB%.db}.${tmp}_${FIX_REPMASK_LAREPEAT_REPEATTRACK}.a2"

        	done > scrub_05_TKcombine_single_${FIX_DB%.db}.${slurmID}.plan
        fi 
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_05_TKcombine_single_${FIX_DB%.db}.${slurmID}.version
    #### TKhomogenize
    elif [[ ${currentStep} -eq 6 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_06_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done             
         # we need the name of the repeat track, especially if the plan starts with step6
        setLArepeatOptions 1
        ### create TKhomogenize commands
        for x in $(seq 0 $((${numRepeatTracks}-1)))
        do 
            tmp=$(echo ${SCRUB_LAREPEAT_OPT[${x}]} | awk '{print $NF}')_${FIX_REPMASK_LAREPEAT_REPEATTRACK}
            for y in $(seq 1 ${fixblocks})
            do 
                echo "${MARVEL_PATH}/bin/TKhomogenize -i  ${tmp} -I h${tmp} -b ${y} ${FIX_DB%.db} ${FIX_DB%.db}.${y}.dalign.las"
            done 
    	done > scrub_06_TKhomogenize_block_${FIX_DB%.db}.${slurmID}.plan
    	echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_06_TKhomogenize_block_${FIX_DB%.db}.${slurmID}.version         
    #### TKcombine
    elif [[ ${currentStep} -eq 7 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_07_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done             
         # we need the name of the repeat track, especially if the plan starts with step7
        setLArepeatOptions 1
        ### find and set TKcombine options
        setTKcombineOptions 0
        ### create TKcombine commands
        for x in $(seq 0 $((${numRepeatTracks}-1)))
        do 
            tmp=$(echo ${SCRUB_LAREPEAT_OPT[${x}]} | awk '{print $NF}')_${FIX_REPMASK_LAREPEAT_REPEATTRACK}
            echo "${MARVEL_PATH}/bin/TKcombine${SCRUB_TKCOMBINE_OPT} ${FIX_DB%.db} h${tmp} \#.h${tmp}"
    	done > scrub_07_TKcombine_single_${FIX_DB%.db}.${slurmID}.plan         
        ### find and set TKcombine options, but ignore the -d flag if it was set
        setTKcombineOptions 1

        echo "${MARVEL_PATH}/bin/TKcombine${SCRUB_TKCOMBINE_OPT} ${FIX_DB%.db} ${FIX_REPMASK_TANMASK_TRACK}_dust ${FIX_REPMASK_TANMASK_TRACK} dust" >> scrub_07_TKcombine_single_${FIX_DB%.db}.${slurmID}.plan
        for x in $(seq 0 $((${numRepeatTracks}-1)))
        do
            tmp=$(echo ${SCRUB_LAREPEAT_OPT[${x}]} | awk '{print $NF}')_${FIX_REPMASK_LAREPEAT_REPEATTRACK}
            echo "${MARVEL_PATH}/bin/TKcombine${SCRUB_TKCOMBINE_OPT} ${FIX_DB%.db} f${tmp} h${tmp} ${tmp}" 
            echo "${MARVEL_PATH}/bin/TKcombine${SCRUB_TKCOMBINE_OPT} ${FIX_DB%.db} f${tmp}_${FIX_REPMASK_TANMASK_TRACK} f${tmp} ${FIX_REPMASK_TANMASK_TRACK}" 
            echo "${MARVEL_PATH}/bin/TKcombine${SCRUB_TKCOMBINE_OPT} ${FIX_DB%.db} f${tmp}_${FIX_REPMASK_TANMASK_TRACK}_dust f${tmp}_${FIX_REPMASK_TANMASK_TRACK} dust"            
		done >> scrub_07_TKcombine_single_${FIX_DB%.db}.${slurmID}.plan
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_07_TKcombine_single_${FIX_DB%.db}.${slurmID}.version      
    #### LAstitch
    elif [[ ${currentStep} -eq 8 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_08_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done             
        ### find and set LAstitch options 
        setLAstitchOptions 1 
        ### create LAstitch commands
        for x in $(seq 1 ${fixblocks})
        do 
            echo "${MARVEL_PATH}/bin/LAstitch${SCRUB_STITCH_OPT} ${FIX_DB%.db} ${FIX_DB%.db}.${x}.dalign.las ${FIX_DB%.db}.${x}.dalignStitch.las"
		done > scrub_08_LAstitch_block_${FIX_DB%.db}.${slurmID}.plan
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_08_LAstitch_block_${FIX_DB%.db}.${slurmID}.version                 
    #### LAq
    elif [[ ${currentStep} -eq 9 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_09_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 

        ### find and set LAq options 
        setLAqOptions
        ### create LAq commands
        for x in $(seq 1 ${fixblocks})
        do 
            echo "${MARVEL_PATH}/bin/LAq${SCRUB_LAQ_OPT} -T trim0_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}_dalign -Q q0_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}_dalign -b ${x} ${FIX_DB%.db} ${FIX_DB%.db}.${x}.dalignStitch.las"
    	done > scrub_09_LAq_block_${FIX_DB%.db}.${slurmID}.plan
    	echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_09_LAq_block_${FIX_DB%.db}.${slurmID}.version               
    #### TKmerge    
    elif [[ ${currentStep} -eq 10 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_10_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done         
        # we need the name of the repeat track, especially if the plan starts with step10
        if [[ -z ${SCRUB_LAQ_OPT} ]]
        then 
            setLAqOptions
        fi
        ### find and set TKmerge options 
        if [[ -z ${SCRUB_TKMERGE_OPT} ]]
        then 
            setTKmergeOptions
        fi
        
        ### create TKmerge command
        echo "${MARVEL_PATH}/bin/TKmerge${SCRUB_TKMERGE_OPT} ${FIX_DB%.db} trim0_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}_dalign" > scrub_10_TKmerge_block_${FIX_DB%.db}.${slurmID}.plan 
        echo "${MARVEL_PATH}/bin/TKmerge${SCRUB_TKMERGE_OPT} ${FIX_DB%.db} q0_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}_dalign" >> scrub_10_TKmerge_block_${FIX_DB%.db}.${slurmID}.plan
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_10_TKmerge_block_${FIX_DB%.db}.${slurmID}.version               
    #### LAgap
    elif [[ ${currentStep} -eq 11 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_11_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done         
        setLAgapOptions 1
        ### create LAq commands
        for x in $(seq 1 ${fixblocks})
        do 
            breakChim=""
            if [[ -n ${FIX_SCRUB_LAGAP_DISCARD_CHIMERS} ]]
            then 
                breakChim=" -C ${FIX_DB%.db}.${x}.dalignGap.chimers.txt"
            fi
            echo "${MARVEL_PATH}/bin/LAgap${SCRUB_LAGAP_OPT} -t trim0_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}_dalign${breakChim} ${FIX_DB%.db} ${FIX_DB%.db}.${x}.dalignStitch.las ${FIX_DB%.db}.${x}.dalignGap.las"
    	done > scrub_11_LAgap_block_${FIX_DB%.db}.${slurmID}.plan      
    	echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_11_LAgap_block_${FIX_DB%.db}.${slurmID}.version         
    #### LAq
    elif [[ ${currentStep} -eq 12 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_12_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 

        ### find and set LAq options 
        setLAqOptions
        ### create LAq commands
        for x in $(seq 1 ${fixblocks})
        do 
            echo "${MARVEL_PATH}/bin/LAq${SCRUB_LAQ_OPT} -u -T trim1_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}_dalign -t trim0_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}_dalign -q q0_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}_dalign -b ${x} ${FIX_DB%.db} ${FIX_DB%.db}.${x}.dalignGap.las"
    done > scrub_12_LAq_block_${FIX_DB%.db}.${slurmID}.plan 
    echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_12_LAq_block_${FIX_DB%.db}.${slurmID}.version              
    #### TKmerge    
    elif [[ ${currentStep} -eq 13 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_13_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done         
        # we need the name of the repeat track, especially if the plan starts with step10
        if [[ -z ${SCRUB_LAQ_OPT} ]]
        then 
            setLAqOptions
        fi
        ### find and set TKmerge options 
        if [[ -z ${SCRUB_TKMERGE_OPT} ]]
        then 
            setTKmergeOptions
        fi
        
        ### create TKmerge command
        echo "${MARVEL_PATH}/bin/TKmerge${SCRUB_TKMERGE_OPT} ${FIX_DB%.db} trim1_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}_dalign" > scrub_13_TKmerge_single_${FIX_DB%.db}.${slurmID}.plan 
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_13_TKmerge_single_${FIX_DB%.db}.${slurmID}.version
	#### LAseparate 
    elif [[ ${currentStep} -eq 14 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_14_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 

        ### find and set repcomp options 
        setLAseparateOptions 0

        for x in $(seq 1 ${fixblocks}); 
        do 
            sdir=$(getSubDirName ${FIX_SCRUB_DALIGNER_RUNID} ${x})
            mkdir -p ${sdir}_ForRepComp
            mkdir -p ${sdir}_NoRepComp
            for y in $(seq 1 ${fixblocks}); 
            do 
                if [[ ! -f ${sdir}/${FIX_DB%.db}.${x}.${FIX_DB%.db}.${y}.las ]]
                then
                    (>&2 echo "step ${currentStep} in FIX_SCRUB_TYPE ${FIX_SCRUB_TYPE}: File missing ${sdir}/${FIX_DB%.db}.${x}.${FIX_DB%.db}.${y}.las!!")
                    exit 1                    
                fi
                echo "${MARVEL_PATH}/bin/LAseparate${SCRUB_LASEPARATE_OPT} ${FIX_DB%.db} ${sdir}/${FIX_DB%.db}.${x}.${FIX_DB%.db}.${y}.las ${sdir}_ForRepComp/${FIX_DB%.db}.${x}.${FIX_DB%.db}.${y}.las ${sdir}_NoRepComp/${FIX_DB%.db}.${x}.${FIX_DB%.db}.${y}.las"                
            done 
		done > scrub_14_LAseparate_block_${FIX_DB%.db}.${slurmID}.plan
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_14_LAseparate_block_${FIX_DB%.db}.${slurmID}.version
    #### repcomp 
    elif [[ ${currentStep} -eq 15 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_15_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        ### find and set repcomp options 
        setRepcompOptions
        cmdLine=1
        for x in $(seq 1 ${fixblocks}); 
        do 
            srcDir=$(getSubDirName ${FIX_SCRUB_DALIGNER_RUNID} ${x})_ForRepComp
            desDir=$(getSubDirName ${FIX_SCRUB_REPCOMP_RUNID} ${x})

            if [[ ! -d ${desDir} ]]
            then
                mkdir -p ${desDir}
            fi
            start=${x}

            for y in $(seq ${start} ${fixblocks}); 
            do 
                movDir=$(getSubDirName ${FIX_SCRUB_REPCOMP_RUNID} ${y})
                if [[ -f ${srcDir}/${FIX_DB%.db}.${x}.${FIX_DB%.db}.${y}.las ]]
                then 
                    if [[ -n ${FIX_SCRUB_REPCOMP_NUMACTL} && ${FIX_SCRUB_REPCOMP_NUMACTL} -gt 0 ]] && [[ "x${SLURM_NUMACTL}" == "x" || ${SLURM_NUMACTL} -eq 0 ]]
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
                    echo -n "${NUMACTL}${REPCOMP_PATH}/bin/repcomp${SCRUB_REPCOMP_OPT} -T/tmp/${FIX_DB%.db}.${x}.${y} ${desDir}/${FIX_DB%.db}.repcomp.${x}.${y} ${FIX_DAZZ_DB%.db} ${srcDir}/${FIX_DB%.db}.${x}.${FIX_DB%.db}.${y}.las"
                    cmdLine=$((${cmdLine}+1))
                    if [[ $x -eq $y ]]
                    then
                        echo ""
                    else    
                        echo " && mv ${desDir}/${FIX_DB%.db}.repcomp.${x}.${y}_r.las ${movDir}/"
                    fi
                else
                    (>&2 echo "step ${currentStep} in FIX_SCRUB_TYPE ${FIX_SCRUB_TYPE}: File missing ${srcDir}/${FIX_DB%.db}.${x}.${FIX_DB%.db}.${y}.las!!")
                    exit 1
                fi
            done 
		done > scrub_15_repcomp_block_${FIX_DB%.db}.${slurmID}.plan
		echo "repcomp $(git --git-dir=${REPCOMP_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_15_repcomp_block_${FIX_DB%.db}.${slurmID}.version  
    #### LAmerge
    elif [[ ${currentStep} -eq 16 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_16_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        
        ### find and set LAmerge options 
        setLAmergeOptions
        setRepcompOptions
        ### create LAmerge commands
        for x in $(seq 1 ${fixblocks})
        do  
            echo "${MARVEL_PATH}/bin/LAmerge${SCRUB_LAMERGE_OPT} ${FIX_DB%.db} ${FIX_DB%.db}.${x}.repcomp.las $(getSubDirName ${FIX_SCRUB_REPCOMP_RUNID} ${x}) $(getSubDirName ${FIX_SCRUB_DALIGNER_RUNID} ${x})_NoRepComp"
		done > scrub_16_LAmerge_block_${FIX_DB%.db}.${slurmID}.plan
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_15_repcomp_block_${FIX_DB%.db}.${slurmID}.version   
    #### LArepeat
    elif [[ ${currentStep} -eq 17 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_17_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 

        setLArepeatOptions 2

        if [[ ${numRepeatTracks} -eq 0 ]]
        then 
            exit 1
        fi    
    
        for x in $(seq 0 $((${numRepeatTracks}-1)))
        do 
            ### create LArepeat commands
            for y in $(seq 1 ${fixblocks})
            do 
                echo "${MARVEL_PATH}/bin/LArepeat${SCRUB_LAREPEAT_OPT[$x]} -b ${y} ${FIX_DB%.db} ${FIX_DB%.db}.${y}.repcomp.las"
            done
		done > scrub_17_LArepeat_block_${FIX_DB%.db}.${slurmID}.plan
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_17_LArepeat_block_${FIX_DB%.db}.${slurmID}.version      
    #### TKmerge 
    elif [[ ${currentStep} -eq 18 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_18_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done         
        # we need the name of the repeat track, especially if the plan starts with step10
        setLArepeatOptions 2
        ### find and set TKmerge options 
        if [[ -z ${SCRUB_TKMERGE_OPT} ]]
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
            echo "${MARVEL_PATH}/bin/TKmerge${SCRUB_TKMERGE_OPT} ${FIX_DB%.db} $(echo ${SCRUB_LAREPEAT_OPT[${x}]} | awk '{print $NF}')" 
		done > scrub_18_TKmerge_single_${FIX_DB%.db}.${slurmID}.plan
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_18_TKmerge_single_${FIX_DB%.db}.${slurmID}.version       
    #### TKcombine 
    elif [[ ${currentStep} -eq 19 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_19_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done             
         # we need the name of the repeat track, especially if the plan starts with step11
        setLArepeatOptions 2
        ### find and set TKcombine options
        setTKcombineOptions 1

        if [[ ${numRepeatTracks} -eq 0 ]]
        then 
            exit 1
        fi 
        ### create TKmerge command
        if [[ -n ${FIX_REPMASK_REPEATTRACK} ]]
        then 
            for x in $(seq 0 $((${numRepeatTracks}-1)))
            do 
                tmp=$(echo ${SCRUB_LAREPEAT_OPT[${x}]} | awk '{print $NF}')
                echo "${MARVEL_PATH}/bin/TKcombine${SCRUB_TKCOMBINE_OPT} ${FIX_DB%.db} ${tmp}_${FIX_REPMASK_LAREPEAT_REPEATTRACK} ${tmp} ${FIX_REPMASK_REPEATTRACK}" 
    		done > scrub_19_TKcombine_single_${FIX_DB%.db}.${slurmID}.plan         
        else
            for x in $(seq 0 $((${numRepeatTracks}-1)))
            do 
                tmp=$(echo ${SCRUB_LAREPEAT_OPT[${x}]} | awk '{print $NF}')
                echo "ln -s -f .${FIX_DB%.db}.${tmp}.d2 .${FIX_DB%.db}.${tmp}_${FIX_REPMASK_LAREPEAT_REPEATTRACK}.d2"
                echo "ln -s -f .${FIX_DB%.db}.${tmp}.a2 .${FIX_DB%.db}.${tmp}_${FIX_REPMASK_LAREPEAT_REPEATTRACK}.a2"

			done > scrub_19_TKcombine_single_${FIX_DB%.db}.${slurmID}.plan
        fi 
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_19_TKcombine_single_${FIX_DB%.db}.${slurmID}.version
    #### TKhomogenize
    elif [[ ${currentStep} -eq 20 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_20_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done             
         # we need the name of the repeat track, especially if the plan starts with step6
        setLArepeatOptions 2
        ### create TKhomogenize commands
        for x in $(seq 0 $((${numRepeatTracks}-1)))
        do 
            tmp=$(echo ${SCRUB_LAREPEAT_OPT[${x}]} | awk '{print $NF}')_${FIX_REPMASK_LAREPEAT_REPEATTRACK}
            for y in $(seq 1 ${fixblocks})
            do 
                echo "${MARVEL_PATH}/bin/TKhomogenize -i  ${tmp} -I h${tmp} -b ${y} ${FIX_DB%.db} ${FIX_DB%.db}.${y}.repcomp.las"
            done 
		done > scrub_20_TKhomogenize_block_${FIX_DB%.db}.${slurmID}.plan
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_20_TKhomogenize_block_${FIX_DB%.db}.${slurmID}.version         
    #### TKcombine
    elif [[ ${currentStep} -eq 21 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_21_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done             
         # we need the name of the repeat track, especially if the plan starts with step7
        setLArepeatOptions 2
        ### find and set TKcombine options
        setTKcombineOptions 0
        ### create TKcombine commands
        for x in $(seq 0 $((${numRepeatTracks}-1)))
        do 
            tmp=$(echo ${SCRUB_LAREPEAT_OPT[${x}]} | awk '{print $NF}')_${FIX_REPMASK_LAREPEAT_REPEATTRACK}
            echo "${MARVEL_PATH}/bin/TKcombine${SCRUB_TKCOMBINE_OPT} ${FIX_DB%.db} h${tmp} \#.h${tmp}"
		done > scrub_21_TKcombine_single_${FIX_DB%.db}.${slurmID}.plan         
        ### find and set TKcombine options, but ignore the -d flag if it was set
        setTKcombineOptions 1

        for x in $(seq 0 $((${numRepeatTracks}-1)))
        do
            tmp=$(echo ${SCRUB_LAREPEAT_OPT[${x}]} | awk '{print $NF}')_${FIX_REPMASK_LAREPEAT_REPEATTRACK}
            echo "${MARVEL_PATH}/bin/TKcombine${SCRUB_TKCOMBINE_OPT} ${FIX_DB%.db} f${tmp} h${tmp} ${tmp}" 
            echo "${MARVEL_PATH}/bin/TKcombine${SCRUB_TKCOMBINE_OPT} ${FIX_DB%.db} f${tmp}_${FIX_REPMASK_TANMASK_TRACK} f${tmp} ${FIX_REPMASK_TANMASK_TRACK}" 
            echo "${MARVEL_PATH}/bin/TKcombine${SCRUB_TKCOMBINE_OPT} ${FIX_DB%.db} f${tmp}_${FIX_REPMASK_TANMASK_TRACK}_dust f${tmp}_${FIX_REPMASK_TANMASK_TRACK} dust"            
		done >> scrub_21_TKcombine_single_${FIX_DB%.db}.${slurmID}.plan
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_21_TKcombine_single_${FIX_DB%.db}.${slurmID}.version        
    #### LAstitch
    elif [[ ${currentStep} -eq 22 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_22_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done             
        ### find and set LAstitch options 
        setLAstitchOptions 2 
        ### create LAstitch commands
        for x in $(seq 1 ${fixblocks})
        do 
            echo "${MARVEL_PATH}/bin/LAstitch${SCRUB_STITCH_OPT} ${FIX_DB%.db} ${FIX_DB%.db}.${x}.repcomp.las ${FIX_DB%.db}.${x}.repcompStitch.las"
		done > scrub_22_LAstitch_block_${FIX_DB%.db}.${slurmID}.plan
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_22_LAstitch_block_${FIX_DB%.db}.${slurmID}.version                 
    #### LAq
    elif [[ ${currentStep} -eq 23 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_23_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 

        ### find and set LAq options 
        setLAqOptions
        ### create LAq commands
        for x in $(seq 1 ${fixblocks})
        do 
            echo "${MARVEL_PATH}/bin/LAq${SCRUB_LAQ_OPT} -T trim0_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}_repcomp -Q q0_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}_repcomp -b ${x} ${FIX_DB%.db} ${FIX_DB%.db}.${x}.repcompStitch.las"
		done > scrub_23_LAq_block_${FIX_DB%.db}.${slurmID}.plan
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_23_LAq_block_${FIX_DB%.db}.${slurmID}.version         
    #### TKmerge    
    elif [[ ${currentStep} -eq 24 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_24_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done         
        # we need the name of the repeat track, especially if the plan starts with step10
        if [[ -z ${SCRUB_LAQ_OPT} ]]
        then 
            setLAqOptions
        fi
        ### find and set TKmerge options 
        if [[ -z ${SCRUB_TKMERGE_OPT} ]]
        then 
            setTKmergeOptions
        fi
        
        ### create TKmerge command
        echo "${MARVEL_PATH}/bin/TKmerge${SCRUB_TKMERGE_OPT} ${FIX_DB%.db} trim0_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}_repcomp" > scrub_24_TKmerge_block_${FIX_DB%.db}.${slurmID}.plan
        echo "${MARVEL_PATH}/bin/TKmerge${SCRUB_TKMERGE_OPT} ${FIX_DB%.db} q0_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}_repcomp" >> scrub_24_TKmerge_block_${FIX_DB%.db}.${slurmID}.plan
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_24_TKmerge_block_${FIX_DB%.db}.${slurmID}.version               
    #### LAgap
    elif [[ ${currentStep} -eq 25 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_25_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done         
        setLAgapOptions 2
        ### create LAq commands
        for x in $(seq 1 ${fixblocks})
        do 
            breakChim=""
            if [[ -n ${FIX_SCRUB_LAGAP_DISCARD_CHIMERS} ]]
            then 
                breakChim=" -C ${FIX_DB%.db}.${x}.repcompGap.chimers.txt"
            fi
            echo "${MARVEL_PATH}/bin/LAgap${SCRUB_LAGAP_OPT} -t trim0_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}_repcomp${breakChim} ${FIX_DB%.db} ${FIX_DB%.db}.${x}.repcompStitch.las ${FIX_DB%.db}.${x}.repcompGap.las"
		done > scrub_25_LAgap_block_${FIX_DB%.db}.${slurmID}.plan
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_25_LAgap_block_${FIX_DB%.db}.${slurmID}.version               
    #### LAq
    elif [[ ${currentStep} -eq 26 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_26_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 

        ### find and set LAq options 
        setLAqOptions
        ### create LAq commands
        for x in $(seq 1 ${fixblocks})
        do 
            echo "${MARVEL_PATH}/bin/LAq${SCRUB_LAQ_OPT} -u -T trim1_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}_repcomp -t trim0_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}_repcomp -q q0_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}_repcomp -b ${x} ${FIX_DB%.db} ${FIX_DB%.db}.${x}.repcompGap.las"
		done > scrub_26_LAq_block_${FIX_DB%.db}.${slurmID}.plan
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_26_LAq_block_${FIX_DB%.db}.${slurmID}.version               
    #### TKmerge    
    elif [[ ${currentStep} -eq 27 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_27_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done         
        # we need the name of the repeat track, especially if the plan starts with step10
        if [[ -z ${SCRUB_LAQ_OPT} ]]
        then 
            setLAqOptions
        fi
        ### find and set TKmerge options 
        if [[ -z ${SCRUB_TKMERGE_OPT} ]]
        then 
            setTKmergeOptions
        fi
        
        ### create TKmerge command
        echo "${MARVEL_PATH}/bin/TKmerge${SCRUB_TKMERGE_OPT} ${FIX_DB%.db} trim1_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}_repcomp" > scrub_27_TKmerge_single_${FIX_DB%.db}.${slurmID}.plan
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_27_TKmerge_single_${FIX_DB%.db}.${slurmID}.version           
    #### LAseparate 
    elif [[ ${currentStep} -eq 28 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_28_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 

        ### find and set repcomp options 
        setRepcompOptions
        setDalignerOptions
        setLAseparateOptions 1
        for x in $(seq 1 ${fixblocks}); 
        do 
            sdir=$(getSubDirName ${FIX_SCRUB_REPCOMP_RUNID} ${x})
            mkdir -p ${sdir}_ForForceAlign
            mkdir -p ${sdir}_NoForceAlign
            for y in $(seq 1 ${fixblocks}); 
            do 
                infile=""
                if [[ $x -le $y ]]
                then    
                    infile=${FIX_DB%.db}.repcomp.${x}.${y}_f.las 

                    if [[ ! -f ${sdir}/${infile} ]]
                    then
                        (>&2 echo "step ${currentStep} in FIX_SCRUB_TYPE ${FIX_SCRUB_TYPE}: File missing ${sdir}/${infile}!!")
                        exit 1                    
                    fi
                    echo "${MARVEL_PATH}/bin/LAseparate${SCRUB_LASEPARATE_OPT} ${FIX_DB%.db} ${sdir}/${infile} ${sdir}_ForForceAlign/${infile} ${sdir}_NoForceAlign/${infile}"                
                fi
                if [[ $x -ge $y ]]
                then    
                    infile=${FIX_DB%.db}.repcomp.${y}.${x}_r.las 

                    if [[ ! -f ${sdir}/${infile} ]]
                    then
                        (>&2 echo "step ${currentStep} in FIX_SCRUB_TYPE ${FIX_SCRUB_TYPE}: File missing ${sdir}/${infile}!!")
                        exit 1                    
                    fi
                    echo "${MARVEL_PATH}/bin/LAseparate${SCRUB_LASEPARATE_OPT} ${FIX_DB%.db} ${sdir}/${infile} ${sdir}_ForForceAlign/${infile} ${sdir}_NoForceAlign/${infile}"                
                fi
            done 
		done > scrub_28_LAseparate_block_${FIX_DB%.db}.${slurmID}.plan   
        for x in $(seq 1 ${fixblocks}); 
        do 
            sdir=$(getSubDirName ${FIX_SCRUB_DALIGNER_RUNID} ${x})_NoRepComp
            mkdir -p ${sdir}_ForForceAlign
            mkdir -p ${sdir}_NoForceAlign
            for y in $(seq 1 ${fixblocks}); 
            do 
                if [[ ! -f ${sdir}/${FIX_DB%.db}.${x}.${FIX_DB%.db}.${y}.las ]]
                then
                    (>&2 echo "step ${currentStep} in FIX_SCRUB_TYPE ${FIX_SCRUB_TYPE}: File missing ${sdir}/${FIX_DB%.db}.${x}.${FIX_DB%.db}.${y}.las!!")
                    exit 1                    
                fi
                echo "${MARVEL_PATH}/bin/LAseparate${SCRUB_LASEPARATE_OPT} ${FIX_DB%.db} ${sdir}/${FIX_DB%.db}.${x}.${FIX_DB%.db}.${y}.las ${sdir}_ForForceAlign/${FIX_DB%.db}.${x}.${FIX_DB%.db}.${y}.las ${sdir}_NoForceAlign/${FIX_DB%.db}.${x}.${FIX_DB%.db}.${y}.las"                
            done             
		done >> scrub_28_LAseparate_block_${FIX_DB%.db}.${slurmID}.plan
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_28_LAseparate_block_${FIX_DB%.db}.${slurmID}.version        
    #### forcealign
    elif [[ ${currentStep} -eq 29 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_29_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        ### find and set forcealign options 
        setForcealignOptions   
        cmdLine=1        
        for x in $(seq 1 ${fixblocks}); 
        do 
            srcDir=$(getSubDirName ${FIX_SCRUB_REPCOMP_RUNID} ${x})_ForForceAlign
            desDir=$(getSubDirName ${FIX_SCRUB_FORCEALIGN_RUNID} ${x})

            if [[ ! -d ${desDir} ]]
            then
                mkdir -p ${desDir}
            fi
            start=${x}

            for y in $(seq ${start} ${fixblocks}); 
            do 
                movDir=$(getSubDirName ${FIX_SCRUB_FORCEALIGN_RUNID} ${y})
                inFile=${srcDir}/${FIX_DB%.db}.repcomp.${x}.${y}_f.las

                if [[ -f ${inFile} ]]
                then 
                    if [[ -n ${FIX_SCRUB_FORCEALIGN_NUMACTL} && ${FIX_SCRUB_FORCEALIGN_NUMACTL} -gt 0 ]] && [[ "x${SLURM_NUMACTL}" == "x" || ${SLURM_NUMACTL} -eq 0 ]]
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
                    echo -n "${NUMACTL}${DACCORD_PATH}/bin/forcealign${SCRUB_FORCEALIGN_OPT} -T/tmp/${FIX_DB%.db}.forcealign.${x}.${y} ${desDir}/${FIX_DB%.db}.forcealign.${x}.${y} ${FIX_DAZZ_DB%.db} ${inFile}"
                    cmdLine=$((${cmdLine}+1))
                    if [[ $x -eq $y ]]
                    then
                        echo ""
                    else    
                        echo " && mv ${desDir}/${FIX_DB%.db}.forcealign.${x}.${y}_r.las ${movDir}"
                    fi
                else
                    (>&2 echo "step ${currentStep} in FIX_SCRUB_TYPE ${FIX_SCRUB_TYPE}: File missing ${inFile}!!")
                    exit 1
                fi
            done 
		done > scrub_29_forcealign_block_${FIX_DB%.db}.${slurmID}.plan
        for x in $(seq 1 ${fixblocks}); 
        do 
            srcDir=$(getSubDirName ${FIX_SCRUB_DALIGNER_RUNID} ${x})_NoRepComp_ForForceAlign
            desDir=$(getSubDirName ${FIX_SCRUB_FORCEALIGN_RUNID} ${x})

            if [[ ! -d ${desDir} ]]
            then
                mkdir -p ${desDir}
            fi
            start=${x}

            for y in $(seq ${start} ${fixblocks}); 
            do 
                movDir=$(getSubDirName ${FIX_SCRUB_FORCEALIGN_RUNID} ${y})
                inFile=${srcDir}/${FIX_DB%.db}.${x}.${FIX_DB%.db}.${y}.las
                
                if [[ -f ${inFile} ]]
                then 
                    if [[ -n ${FIX_SCRUB_FORCEALIGN_NUMACTL} && ${FIX_SCRUB_FORCEALIGN_NUMACTL} -gt 0 ]] && [[ "x${SLURM_NUMACTL}" == "x" || ${SLURM_NUMACTL} -eq 0 ]]
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
                    echo -n "${NUMACTL}${DACCORD_PATH}/bin/forcealign${SCRUB_FORCEALIGN_OPT} -T/tmp/${FIX_DB%.db}.norepcomp.forcealign.${x}.${y} ${desDir}/${FIX_DB%.db}.norepcomp.forcealign.${x}.${y} ${FIX_DAZZ_DB%.db} ${inFile}"
                    cmdLine=$((${cmdLine}+1))
                    if [[ $x -eq $y ]]
                    then
                        echo ""
                    else    
                        echo " && mv ${desDir}/${FIX_DB%.db}.norepcomp.forcealign.${x}.${y}_r.las ${movDir}"
                    fi
                else
                    (>&2 echo "step ${currentStep} in FIX_SCRUB_TYPE ${FIX_SCRUB_TYPE}: File missing ${inFile}!!")
                    exit 1
                fi
            done 
		done >> scrub_29_forcealign_block_${FIX_DB%.db}.${slurmID}.plan
		echo "forcealign $(git --git-dir=${DACCORD_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_29_forcealign_block_${FIX_DB%.db}.${slurmID}.version
    #### LAmerge
    elif [[ ${currentStep} -eq 30 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_30_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        
        ### find and set LAmerge options 
        setLAmergeOptions
        setRepcompOptions
        setForcealignOptions
        ### create LAmerge commands
        for x in $(seq 1 ${fixblocks})
        do  
            echo "${MARVEL_PATH}/bin/LAmerge${SCRUB_LAMERGE_OPT} ${FIX_DB%.db} ${FIX_DB%.db}.${x}.forcealign.las $(getSubDirName ${FIX_SCRUB_FORCEALIGN_RUNID} ${x}) $(getSubDirName ${FIX_SCRUB_REPCOMP_RUNID} ${x})_NoForceAlign $(getSubDirName ${FIX_SCRUB_DALIGNER_RUNID} ${x})_NoRepComp_NoForceAlign"
		done > scrub_30_LAmerge_block_${FIX_DB%.db}.${slurmID}.plan
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_30_LAmerge_block_${FIX_DB%.db}.${slurmID}.version   
    #### LArepeat
    elif [[ ${currentStep} -eq 31 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_31_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 

        setLArepeatOptions 3
        if [[ ${numRepeatTracks} -eq 0 ]]
        then 
            exit 1
        fi    
    
        for x in $(seq 0 $((${numRepeatTracks}-1)))
        do 
            ### create LArepeat commands
            for y in $(seq 1 ${fixblocks})
            do 
                echo "${MARVEL_PATH}/bin/LArepeat${SCRUB_LAREPEAT_OPT[$x]} -b ${y} ${FIX_DB%.db} ${FIX_DB%.db}.${y}.forcealign.las"
            done
		done > scrub_31_LArepeat_block_${FIX_DB%.db}.${slurmID}.plan
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_31_LArepeat_block_${FIX_DB%.db}.${slurmID}.version      
    #### TKmerge 
    elif [[ ${currentStep} -eq 32 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_32_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done         
        # we need the name of the repeat track, especially if the plan starts with step10
        setLArepeatOptions 3
        ### find and set TKmerge options 
        if [[ -z ${SCRUB_TKMERGE_OPT} ]]
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
            echo "${MARVEL_PATH}/bin/TKmerge${SCRUB_TKMERGE_OPT} ${FIX_DB%.db} $(echo ${SCRUB_LAREPEAT_OPT[${x}]} | awk '{print $NF}')" 
		done > scrub_32_TKmerge_single_${FIX_DB%.db}.${slurmID}.plan
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_32_TKmerge_single_${FIX_DB%.db}.${slurmID}.version       
    #### TKcombine 
    elif [[ ${currentStep} -eq 33 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_33_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done             
         # we need the name of the repeat track, especially if the plan starts with step11
        setLArepeatOptions 3
        ### find and set TKcombine options
        setTKcombineOptions 1

        if [[ ${numRepeatTracks} -eq 0 ]]
        then 
            exit 1
        fi 
        ### create TKmerge command
        if [[ -n ${FIX_REPMASK_REPEATTRACK} ]]
        then 
            for x in $(seq 0 $((${numRepeatTracks}-1)))
            do 
                tmp=$(echo ${SCRUB_LAREPEAT_OPT[${x}]} | awk '{print $NF}')
                echo "${MARVEL_PATH}/bin/TKcombine${SCRUB_TKCOMBINE_OPT} ${FIX_DB%.db} ${tmp}_${FIX_REPMASK_LAREPEAT_REPEATTRACK} ${tmp} ${FIX_REPMASK_REPEATTRACK}" 
    		done > scrub_33_TKcombine_single_${FIX_DB%.db}.${slurmID}.plan         
        else
            for x in $(seq 0 $((${numRepeatTracks}-1)))
            do 
                tmp=$(echo ${SCRUB_LAREPEAT_OPT[${x}]} | awk '{print $NF}')
                echo "ln -s -f .${FIX_DB%.db}.${tmp}.d2 .${FIX_DB%.db}.${tmp}_${FIX_REPMASK_LAREPEAT_REPEATTRACK}.d2"
                echo "ln -s -f .${FIX_DB%.db}.${tmp}.a2 .${FIX_DB%.db}.${tmp}_${FIX_REPMASK_LAREPEAT_REPEATTRACK}.a2"

    		done > scrub_33_TKcombine_single_${FIX_DB%.db}.${slurmID}.plan
        fi
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_33_TKcombine_single_${FIX_DB%.db}.${slurmID}.version 
    #### TKhomogenize
    elif [[ ${currentStep} -eq 34 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_34_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done             
         # we need the name of the repeat track, especially if the plan starts with step6
        setLArepeatOptions 3
        ### create TKhomogenize commands
        for x in $(seq 0 $((${numRepeatTracks}-1)))
        do 
            tmp=$(echo ${SCRUB_LAREPEAT_OPT[${x}]} | awk '{print $NF}')_${FIX_REPMASK_LAREPEAT_REPEATTRACK}
            for y in $(seq 1 ${fixblocks})
            do 
                echo "${MARVEL_PATH}/bin/TKhomogenize -i  ${tmp} -I h${tmp} -b ${y} ${FIX_DB%.db} ${FIX_DB%.db}.${y}.forcealign.las"
            done 
		done > scrub_34_TKhomogenize_block_${FIX_DB%.db}.${slurmID}.plan
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_34_TKhomogenize_block_${FIX_DB%.db}.${slurmID}.version         
    #### TKcombine
    elif [[ ${currentStep} -eq 35 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_35_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done             
         # we need the name of the repeat track, especially if the plan starts with step7
        setLArepeatOptions 3
        ### find and set TKcombine options
        setTKcombineOptions 0
        ### create TKcombine commands
        for x in $(seq 0 $((${numRepeatTracks}-1)))
        do 
            tmp=$(echo ${SCRUB_LAREPEAT_OPT[${x}]} | awk '{print $NF}')_${FIX_REPMASK_LAREPEAT_REPEATTRACK}
            echo "${MARVEL_PATH}/bin/TKcombine${SCRUB_TKCOMBINE_OPT} ${FIX_DB%.db} h${tmp} \#.h${tmp}"
		done > scrub_35_TKcombine_single_${FIX_DB%.db}.${slurmID}.plan         
        ### find and set TKcombine options, but ignore the -d flag if it was set
        setTKcombineOptions 1

        for x in $(seq 0 $((${numRepeatTracks}-1)))
        do
            tmp=$(echo ${SCRUB_LAREPEAT_OPT[${x}]} | awk '{print $NF}')_${FIX_REPMASK_LAREPEAT_REPEATTRACK}
            echo "${MARVEL_PATH}/bin/TKcombine${SCRUB_TKCOMBINE_OPT} ${FIX_DB%.db} f${tmp} h${tmp} ${tmp}" 
            echo "${MARVEL_PATH}/bin/TKcombine${SCRUB_TKCOMBINE_OPT} ${FIX_DB%.db} f${tmp}_${FIX_REPMASK_TANMASK_TRACK} f${tmp} ${FIX_REPMASK_TANMASK_TRACK}" 
            echo "${MARVEL_PATH}/bin/TKcombine${SCRUB_TKCOMBINE_OPT} ${FIX_DB%.db} f${tmp}_${FIX_REPMASK_TANMASK_TRACK}_dust f${tmp}_${FIX_REPMASK_TANMASK_TRACK} dust"            
		done >> scrub_35_TKcombine_single_${FIX_DB%.db}.${slurmID}.plan
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_35_TKcombine_single_${FIX_DB%.db}.${slurmID}.version        
    #### LAstitch
    elif [[ ${currentStep} -eq 36 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_36_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done             
        ### find and set LAstitch options 
        setLAstitchOptions 3 
        ### create LAstitch commands
        for x in $(seq 1 ${fixblocks})
        do 
            echo "${MARVEL_PATH}/bin/LAstitch${SCRUB_STITCH_OPT} ${FIX_DB%.db} ${FIX_DB%.db}.${x}.forcealign.las ${FIX_DB%.db}.${x}.forcealignStitch.las"
		done > scrub_36_LAstitch_block_${FIX_DB%.db}.${slurmID}.plan
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_36_LAstitch_block_${FIX_DB%.db}.${slurmID}.version                 
    #### LAq
    elif [[ ${currentStep} -eq 37 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_37_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 

        ### find and set LAq options 
        setLAqOptions
        ### create LAq commands
        for x in $(seq 1 ${fixblocks})
        do 
            echo "${MARVEL_PATH}/bin/LAq${SCRUB_LAQ_OPT} -T trim0_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}_forcealign -Q q0_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}_forcealign -b ${x} ${FIX_DB%.db} ${FIX_DB%.db}.${x}.forcealignStitch.las"
		done > scrub_37_LAq_block_${FIX_DB%.db}.${slurmID}.plan
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_37_LAq_block_${FIX_DB%.db}.${slurmID}.version               
    #### TKmerge    
    elif [[ ${currentStep} -eq 38 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_38_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done         
        # we need the name of the repeat track, especially if the plan starts with step10
        if [[ -z ${SCRUB_LAQ_OPT} ]]
        then 
            setLAqOptions
        fi
        ### find and set TKmerge options 
        if [[ -z ${SCRUB_TKMERGE_OPT} ]]
        then 
            setTKmergeOptions
        fi
        
        ### create TKmerge command
        echo "${MARVEL_PATH}/bin/TKmerge${SCRUB_TKMERGE_OPT} ${FIX_DB%.db} trim0_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}_forcealign" > scrub_38_TKmerge_block_${FIX_DB%.db}.${slurmID}.plan
        echo "${MARVEL_PATH}/bin/TKmerge${SCRUB_TKMERGE_OPT} ${FIX_DB%.db} q0_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}_forcealign" >> scrub_38_TKmerge_block_${FIX_DB%.db}.${slurmID}.plan       
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_38_TKmerge_block_${FIX_DB%.db}.${slurmID}.version               
    #### LAgap
    elif [[ ${currentStep} -eq 39 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_39_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done         
        setLAgapOptions 3
        ### create LAq commands
        for x in $(seq 1 ${fixblocks})
        do 
            breakChim=""
            if [[ -n ${FIX_SCRUB_LAGAP_DISCARD_CHIMERS} ]]
            then 
                breakChim=" -C ${FIX_DB%.db}.${x}.forcealignGap.chimers.txt"
            fi
            echo "${MARVEL_PATH}/bin/LAgap${SCRUB_LAGAP_OPT} -t trim0_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}_forcealign${breakChim} ${FIX_DB%.db} ${FIX_DB%.db}.${x}.forcealignStitch.las ${FIX_DB%.db}.${x}.forcealignGap.las"
		done > scrub_39_LAgap_block_${FIX_DB%.db}.${slurmID}.plan
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_39_LAgap_block_${FIX_DB%.db}.${slurmID}.version               
    #### LAq
    elif [[ ${currentStep} -eq 40 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_40_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 

        ### find and set LAq options 
        setLAqOptions
        ### create LAq commands
        for x in $(seq 1 ${fixblocks})
        do 
            echo "${MARVEL_PATH}/bin/LAq${SCRUB_LAQ_OPT} -u -T trim1_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}_forcealign -t trim0_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}_forcealign -q q0_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}_forcealign -b ${x} ${FIX_DB%.db} ${FIX_DB%.db}.${x}.forcealignGap.las"
		done > scrub_40_LAq_block_${FIX_DB%.db}.${slurmID}.plan
		echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_40_LAq_block_${FIX_DB%.db}.${slurmID}.version               
    #### TKmerge    
    elif [[ ${currentStep} -eq 41 ]]
    then
        ### clean up plans 
        for x in $(ls scrub_41_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done         
        # we need the name of the repeat track, especially if the plan starts with step10
        if [[ -z ${SCRUB_LAQ_OPT} ]]
        then 
            setLAqOptions
        fi
        ### find and set TKmerge options 
        if [[ -z ${SCRUB_TKMERGE_OPT} ]]
        then 
            setTKmergeOptions
        fi
        
        ### create TKmerge command
        echo "${MARVEL_PATH}/bin/TKmerge${SCRUB_TKMERGE_OPT} ${FIX_DB%.db} trim1_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}_forcealign" > scrub_41_TKmerge_single_${FIX_DB%.db}.${slurmID}.plan
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > scrub_41_TKmerge_single_${FIX_DB%.db}.${slurmID}.version           
    else
        (>&2 echo "step ${currentStep} in FIX_SCRUB_TYPE ${FIX_SCRUB_TYPE} not supported")
        (>&2 echo "valid steps are: ${myTypes[${FIX_SCRUB_TYPE}]}")
        exit 1            
    fi
else
	(>&2echo "unknown FIX_SCRUB_TYPE ${FIX_SCRUB_TYPE}")    
    (>&2 echo "supported types")
    x=0; while [ $x -lt ${#myTypes[*]} ]; do (>&2 echo "${myTypes[${x}]}"); done
    exit 1
fi

exit 0
