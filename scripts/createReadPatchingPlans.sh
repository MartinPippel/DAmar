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

if [[ ! -n "${RAW_PATCH_TYPE}" ]]
then 
    (>&2 echo "cannot create read patching scripts if variable RAW_PATCH_TYPE is not set.")
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

if [[ ${RAW_PATCH_TYPE} -gt 1 && ! -f ${RAW_DAZZ_DB%.db}.db ]]
then 
    (>&2 echo "raw dazzler database ${RAW_DB%.db}.db missing")
    exit 1
fi

function getNumOfDbBlocks()
{
    if [[ ! -f ${RAW_DB%.db}.db ]]
    then
        (>&2 echo "raw database ${RAW_DB%.db}.db not found")
        exit 1
    fi

    blocks=$(grep block ${RAW_DB%.db}.db | awk '{print $3}')
    if [[ ! -n $blocks ]]
    then 
        (>&2 echo "raw database ${RAW_DB%.db}.db has not been partitioned. Run DBsplit first!")
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

function setDalignerOptions()
{
    FIX_DALIGNER_OPT=""
    if [[ -n ${RAW_FIX_DALIGNER_IDENTITY_OVLS} && ${RAW_FIX_DALIGNER_IDENTITY_OVLS} -gt 0 ]]
    then
        FIX_DALIGNER_OPT="${FIX_DALIGNER_OPT} -I"
    fi
    if [[ -n ${RAW_FIX_DALIGNER_KMER} && ${RAW_FIX_DALIGNER_KMER} -gt 0 ]]
    then
        FIX_DALIGNER_OPT="${FIX_DALIGNER_OPT} -k${RAW_FIX_DALIGNER_KMER}"
    fi
    if [[ -n ${RAW_FIX_DALIGNER_ERR} ]]
    then
        FIX_DALIGNER_OPT="${FIX_DALIGNER_OPT} -e${RAW_FIX_DALIGNER_ERR}"
    fi
    if [[ -n ${RAW_FIX_DALIGNER_BIAS} && ${RAW_FIX_DALIGNER_BIAS} -eq 1 ]]
    then
        FIX_DALIGNER_OPT="${FIX_DALIGNER_OPT} -b"
    fi
    if [[ -n ${RAW_FIX_DALIGNER_VERBOSE} && ${RAW_FIX_DALIGNER_VERBOSE} -ne 0 ]]
    then
        FIX_DALIGNER_OPT="${FIX_DALIGNER_OPT} -v"
    fi
    if [[ -n ${RAW_FIX_DALIGNER_TRACESPACE} && ${RAW_FIX_DALIGNER_TRACESPACE} -gt 0 ]]
    then
        FIX_DALIGNER_OPT="${FIX_DALIGNER_OPT} -s${RAW_FIX_DALIGNER_TRACESPACE}"
    fi
    if [[ -n ${RAW_FIX_DALIGNER_ASYMMETRIC} && ${RAW_FIX_DALIGNER_ASYMMETRIC} -ne 0 ]]
    then
        FIX_DALIGNER_OPT="${FIX_DALIGNER_OPT} -A"
    fi    
    if [[ -n ${RAW_FIX_DALIGNER_T} ]]
    then
        FIX_DALIGNER_OPT="${FIX_DALIGNER_OPT} -t${RAW_FIX_DALIGNER_T}"
    fi  
    if [[ -n ${RAW_FIX_DALIGNER_MASK} ]]
    then
        for x in ${RAW_FIX_DALIGNER_MASK}
        do 
            FIX_DALIGNER_OPT="${FIX_DALIGNER_OPT} -m${x}"
        done
    fi
    if [[ -n ${THREADS_daligner} ]]
    then 
        FIX_DALIGNER_OPT="${FIX_DALIGNER_OPT} -T${THREADS_daligner}"
    fi
    if [ ! -n ${RAW_FIX_DALIGNER_DAL} ]
    then
        RAW_FIX_DALIGNER_DAL=8
    fi 
}

function setRepcompOptions()
{
    FIX_REPCOMP_OPT=""
    if [[ -n ${RAW_FIX_REPCOMP_TRACESPACE} && ${RAW_FIX_REPCOMP_TRACESPACE} -gt 0 ]]
    then 
        FIX_REPCOMP_OPT="${FIX_REPCOMP_OPT} --tspace${RAW_FIX_REPCOMP_TRACESPACE}"
    fi
    if [[ -n ${RAW_FIX_REPCOMP_INBLOCKSIZE} ]]
    then 
        FIX_REPCOMP_OPT="${FIX_REPCOMP_OPT} -i${RAW_FIX_REPCOMP_INBLOCKSIZE}"
    fi
    if [[ -n ${RAW_FIX_REPCOMP_KMER} && ${RAW_FIX_REPCOMP_KMER} -gt 0 ]]
    then 
        FIX_REPCOMP_OPT="${FIX_REPCOMP_OPT} -k${RAW_FIX_REPCOMP_KMER}"
    fi
    if [[ -n ${RAW_FIX_REPCOMP_MEM} ]]
    then 
        FIX_REPCOMP_OPT="${FIX_REPCOMP_OPT} -M${RAW_FIX_REPCOMP_MEM}"
    fi
    if [[ -n ${RAW_FIX_REPCOMP_THREADS} && ${RAW_FIX_REPCOMP_THREADS} -gt 0 ]]
    then 
        FIX_REPCOMP_OPT="${FIX_REPCOMP_OPT} -t${RAW_FIX_REPCOMP_THREADS}"
    fi 
    if [[ -n ${RAW_FIX_REPCOMP_CORRELATION} ]]
    then 
        FIX_REPCOMP_OPT="${FIX_REPCOMP_OPT} -e${RAW_FIX_REPCOMP_CORRELATION}"
    fi 
    if [[ -n ${RAW_FIX_REPCOMP_MASK} ]]
    then
        for x in ${RAW_FIX_REPCOMP_MASK}
        do 
            FIX_REPCOMP_OPT="${FIX_REPCOMP_OPT} -m${x}"
        done
    fi
        if [[ -n ${RAW_FIX_REPCOMP_OLEN} && ${RAW_FIX_REPCOMP_OLEN} -gt 0 ]]
    then 
        FIX_REPCOMP_OPT="${FIX_REPCOMP_OPT} -l${RAW_FIX_REPCOMP_OLEN}"
    fi 


    if [[ -z ${RAW_FIX_REPCOMP_RUNID} || ${RAW_FIX_REPCOMP_RUNID} == ${RAW_FIX_DALIGNER_RUNID} ]]
    then
        RAW_FIX_REPCOMP_RUNID=$((${RAW_FIX_DALIGNER_RUNID}+1))
    fi
}

function setLAmergeOptions()
{
    FIX_LAMERGE_OPT=""
    if [[ -n ${RAW_FIX_LAMERGE_NFILES} && ${RAW_FIX_LAMERGE_NFILES} -gt 0 ]]
    then
        FIX_LAMERGE_OPT="${FIX_LAMERGE_OPT} -n ${RAW_FIX_LAMERGE_NFILES}"
    fi
}

function setLArepeatOptions()
{
    ### find and set LArepeat options 
    FIX_LAREPEAT_OPT=""
    if [[ -n ${RAW_FIX_LAREPEAT_LEAVE_COV} ]]
    then
        FIX_LAREPEAT_OPT="${FIX_LAREPEAT_OPT} -l ${RAW_FIX_LAREPEAT_LEAVE_COV}"
    else  # set default value
        RAW_FIX_LAREPEAT_LEAVE_COV=1.7
    fi
    if [[ -n ${RAW_FIX_LAREPEAT_ENTER_COV} ]]
    then
        FIX_LAREPEAT_OPT="${FIX_LAREPEAT_OPT} -h ${RAW_FIX_LAREPEAT_ENTER_COV}"
    else  # set default value
        RAW_FIX_LAREPEAT_ENTER_COV=2.0
    fi
    if [[ -n ${RAW_FIX_LAREPEAT_COV} ]]
    then
        FIX_LAREPEAT_OPT="${FIX_LAREPEAT_OPT} -c ${RAW_FIX_LAREPEAT_COV}"
    fi
    if [[ -n ${RAW_FIX_LAREPEAT_REPEATTRACK} ]]
    then
        FIX_LAREPEAT_OPT="${FIX_LAREPEAT_OPT} -t ${RAW_FIX_LAREPEAT_REPEATTRACK}"
    else
        ptype=""
        if [[ $1 -eq 1 ]]
        then 
            ptype="_dalign"
        elif [[ $1 -eq 2 ]]
        then 
            ptype="_repcomp"
        elif [[ $1 -eq 3 ]]
        then 
            ptype="_forcealign"        
        else
            (>&2 echo "Unknown RAW_PATCH_TYPE=${RAW_PATCH_TYPE} !!!")
            exit 1            
        fi 
        
        RAW_DAZZ_FIX_LAREPEAT_THRESHOLD=$(echo "${RAW_COV} ${RAW_FIX_LAREPEAT_ENTER_COV}" | awk '{printf "%d", $1*$2}')
        RAW_DAZZ_FIX_LAREPEAT_REPEATTRACK=repeats_c${RAW_DAZZ_FIX_LAREPEAT_THRESHOLD}${ptype}      

        if [[ -n ${RAW_FIX_LAREPEAT_COV} ]]
        then
            RAW_FIX_LAREPEAT_REPEATTRACK=repeats_c${RAW_FIX_LAREPEAT_COV}_l${RAW_FIX_LAREPEAT_LEAVE_COV}h${RAW_FIX_LAREPEAT_ENTER_COV}${ptype}
        else
        	if [[ -n ${RAW_FIX_LAREPEAT_MAX_COV} && ${RAW_FIX_LAREPEAT_MAX_COV} -gt 100 ]]
        	then 
        		FIX_LAREPEAT_OPT="${FIX_LAREPEAT_OPT} -M ${RAW_FIX_LAREPEAT_MAX_COV}"
			elif [[ -n ${RAW_COV} && $((${RAW_COV}+20)) -gt 100 ]]
			then
				FIX_LAREPEAT_OPT="${FIX_LAREPEAT_OPT} -M 200"
        	fi
        	
            RAW_FIX_LAREPEAT_REPEATTRACK=repeats_calCov_l${RAW_FIX_LAREPEAT_LEAVE_COV}h${RAW_FIX_LAREPEAT_ENTER_COV}${ptype}
        fi
        FIX_LAREPEAT_OPT="${FIX_LAREPEAT_OPT} -t ${RAW_FIX_LAREPEAT_REPEATTRACK}"
    fi
}

function setTKmergeOptions() 
{
    FIX_TKMERGE_OPT=""
    if [[ -n ${RAW_FIX_TKMERGE_DELETE} && ${RAW_FIX_TKMERGE_DELETE} -ne 0 ]]
    then
        FIX_TKMERGE_OPT="${FIX_TKMERGE_OPT} -d"
    fi
}

function setTKcombineOptions() 
{
    ignoreDelete=$1
    FIX_TKCOMBINE_OPT=""
    if [[ ${ignoreDelete} -eq 0 && -n ${RAW_FIX_TKCOMBINE_DELETE} && ${RAW_FIX_TKCOMBINE_DELETE} -ne 0 ]]
    then
        FIX_TKCOMBINE_OPT="${FIX_TKCOMBINE_OPT} -d"
    fi
    if [[ -n ${RAW_FIX_TKCOMBINE_VERBOSE} && ${RAW_FIX_TKCOMBINE_VERBOSE} -ne 0 ]]
    then
        FIX_TKCOMBINE_OPT="${FIX_TKCOMBINE_OPT} -v"
    fi
}

function setLAqOptions()
{
    FIX_LAQ_OPT=""
    adaptQTRIMCUTOFF=""    

    if [[ -n ${RAW_FIX_LAQ_MINSEG} && ${RAW_FIX_LAQ_MINSEG} -ne 0 ]]
    then
        FIX_LAQ_OPT="${FIX_LAQ_OPT} -s ${RAW_FIX_LAQ_MINSEG}"
    else 
        RAW_FIX_LAQ_MINSEG=25
        FIX_LAQ_OPT="${FIX_LAQ_OPT} -s ${RAW_FIX_LAQ_MINSEG}"
    fi

    if [[ -n ${RAW_FIX_LAQ_QTRIMCUTOFF} && ${RAW_FIX_LAQ_QTRIMCUTOFF} -ne 0 ]]
    then
        if [[ -n ${RAW_FIX_DALIGNER_TRACESPACE} && ${RAW_FIX_DALIGNER_TRACESPACE} -ne 100 ]]
        then 
            adaptQTRIMCUTOFF=$(echo "${RAW_FIX_LAQ_QTRIMCUTOFF}*${RAW_FIX_DALIGNER_TRACESPACE}/100+1" | bc)
            FIX_LAQ_OPT="${FIX_LAQ_OPT} -d ${adaptQTRIMCUTOFF}"
        else
            adaptQTRIMCUTOFF=${RAW_FIX_LAQ_QTRIMCUTOFF}
            FIX_LAQ_OPT="${FIX_LAQ_OPT} -d ${adaptQTRIMCUTOFF}"            
        fi
    else 
        if [[ -n ${RAW_FIX_DALIGNER_TRACESPACE} && ${RAW_FIX_DALIGNER_TRACESPACE} -ne 100 ]]
        then 
            RAW_FIX_LAQ_QTRIMCUTOFF=25
            adaptQTRIMCUTOFF=$(echo "${RAW_FIX_LAQ_QTRIMCUTOFF}*${RAW_FIX_DALIGNER_TRACESPACE}/100+1" | bc)
            FIX_LAQ_OPT="${FIX_LAQ_OPT} -d ${adaptQTRIMCUTOFF}"
        else
            adaptQTRIMCUTOFF=25
            RAW_FIX_LAQ_QTRIMCUTOFF=25
            FIX_LAQ_OPT="${FIX_LAQ_OPT} -d ${adaptQTRIMCUTOFF}"            
        fi
    fi
}

## 1st argument ptype
function setLAfixOptions()
{
	ptype=$1
	
	if [[ -z ${RAW_FIX_LAFIX_PATH} ]]
	then
		RAW_FIX_LAFIX_PATH=patchedReads_${ptype}
	else
		RAW_FIX_LAFIX_PATH=${RAW_FIX_LAFIX_PATH}_${ptype}
	fi

	FIX_LAFIX_OPT=""
    if [[ -n ${RAW_FIX_LAFIX_GAP} && ${RAW_FIX_LAQ_MINSEG} -ne 0 ]]
    then
        FIX_LAFIX_OPT="${FIX_LAFIX_OPT} -g ${RAW_FIX_LAFIX_GAP}"
    fi
    if [[ -n ${RAW_FIX_LAFIX_MLEN} && ${RAW_FIX_LAFIX_MLEN} -ne 0 ]]
    then
        FIX_LAFIX_OPT="${FIX_LAFIX_OPT} -x ${RAW_FIX_LAFIX_MLEN}"
    fi
    if [[ -n ${RAW_FIX_LAFIX_LOW_COVERAGE} && ${RAW_FIX_LAFIX_LOW_COVERAGE} -ne 0 ]]
    then
        FIX_LAFIX_OPT="${FIX_LAFIX_OPT} -l"
    fi
    if [[ -n ${RAW_FIX_LAFIX_USEREPEAT} ]]
    then
    	FIX_LAFIX_OPT="${FIX_LAFIX_OPT} -rrepeats_c${RAW_COV}_l${RAW_FIX_LAREPEAT_LEAVE_COV}h${RAW_FIX_LAREPEAT_ENTER_COV}_${ptype}_${RAW_REPMASK_LAREPEAT_REPEATTRACK}_${RAW_REPMASK_TANMASK_TRACK}_dust"        
    fi
    if [[ -n ${RAW_FIX_LAFIX_MAXCHIMERLEN} && ${RAW_FIX_LAFIX_MAXCHIMERLEN} -ne 0 ]]
    then
        FIX_LAFIX_OPT="${FIX_LAFIX_OPT} -C${RAW_FIX_LAFIX_MAXCHIMERLEN}"
    fi
    if [[ -n ${RAW_FIX_LAFIX_MINCHIMERBORDERCOV} && ${RAW_FIX_LAFIX_MINCHIMERBORDERCOV} -ne 0 ]]
    then
        FIX_LAFIX_OPT="${FIX_LAFIX_OPT} -b${RAW_FIX_LAFIX_MINCHIMERBORDERCOV}"
    fi

    if [[ -z ${FIX_LAQ_OPT} ]]
    then
        setLAqOptions
    fi
    if [[ -n ${RAW_FIX_LAFIX_AGGCHIMERDETECT} && ${RAW_FIX_LAFIX_AGGCHIMERDETECT} -ne 0 ]]
    then
        FIX_LAFIX_OPT="${FIX_LAFIX_OPT} -a"
    fi
    if [[ -n ${RAW_FIX_LAFIX_DISCARDCHIMERS} && ${RAW_FIX_LAFIX_DISCARDCHIMERS} -ne 0 ]]
    then
        FIX_LAFIX_OPT="${FIX_LAFIX_OPT} -d"
    fi
    
    FIX_LAFIX_OPT="${FIX_LAFIX_OPT} -q q0_d${RAW_FIX_LAQ_QTRIMCUTOFF}_s${RAW_FIX_LAQ_MINSEG}_${ptype}"

    if [[ -n ${RAW_FIX_LAFIX_TRIM} && ${RAW_FIX_LAFIX_TRIM} -ne 0 ]]
    then
        FIX_LAFIX_OPT="${FIX_LAFIX_OPT} -t trim0_d${RAW_FIX_LAQ_QTRIMCUTOFF}_s${RAW_FIX_LAQ_MINSEG}_${ptype}"
    fi
    
    if [[ -n ${RAW_FIX_LAFIX_FIXCHIMERS} && ${RAW_FIX_LAFIX_FIXCHIMERS} -ne 0 ]]
    then
        FIX_LAFIX_OPT="${FIX_LAFIX_OPT} -X"
    fi
    
    if [[ -n ${RAW_FIX_LAFIX_CONVERTRACKS} ]]
    then
        for x in ${RAW_FIX_LAFIX_CONVERTRACKS}
        do
            FIX_LAFIX_OPT="${FIX_LAFIX_OPT} -c $x"
        done
    fi
}

function setForcealignOptions()
{
    FIX_FORCEALIGN_OPT=""
    if [[ -n ${RAW_FIX_FORCEALIGN_PARTIAL} && ${RAW_FIX_FORCEALIGN_PARTIAL} -ne 0 ]]
    then
        FIX_FORCEALIGN_OPT="${FIX_FORCEALIGN_OPT} --partial"
    fi
    if [[ -n ${RAW_FIX_FORCEALIGN_THREADS} && ${RAW_FIX_FORCEALIGN_THREADS} -gt 0 ]]
    then 
        FIX_FORCEALIGN_OPT="${FIX_FORCEALIGN_OPT} -t${RAW_FIX_FORCEALIGN_THREADS}"
    fi 
    if [[ -n ${RAW_FIX_FORCEALIGN_MAXDIST} && ${RAW_FIX_FORCEALIGN_MAXDIST} -gt 0 ]]
    then 
        FIX_FORCEALIGN_OPT="${FIX_FORCEALIGN_OPT} --maxdist${RAW_FIX_FORCEALIGN_MAXDIST}"
    fi 
    if [[ -n ${RAW_FIX_FORCEALIGN_BORDER} && ${RAW_FIX_FORCEALIGN_BORDER} -gt 0 ]]
    then 
        FIX_FORCEALIGN_OPT="${FIX_FORCEALIGN_OPT} --border${RAW_FIX_FORCEALIGN_BORDER}"
    fi 
    if [[ -n ${RAW_FIX_FORCEALIGN_CORRELATION} ]]
    then 
        FIX_FORCEALIGN_OPT="${FIX_FORCEALIGN_OPT} --correlation${RAW_FIX_FORCEALIGN_CORRELATION}"
    fi 


    if [[ -z ${RAW_FIX_FORCEALIGN_RUNID} ]]
    then
        if [[ -z ${RAW_FIX_REPCOMP_RUNID} ]]
        then
            RAW_FIX_FORCEALIGN_RUNID=$((${RAW_FIX_DALIGNER_RUNID}+2))
        else 
            RAW_FIX_FORCEALIGN_RUNID=$((${RAW_FIX_REPCOMP_RUNID}+1))
        fi
    fi
}

# first argument LAseparate type
function setLAseparateOptions()
{
    FIX_LASEPARATE_OPT=""
    if [[ -n ${RAW_FIX_LASEPARATE_OLEN} && ${RAW_FIX_LASEPARATE_OLEN} -gt 0 ]]
    then 
        FIX_LASEPARATE_OPT="${FIX_LASEPARATE_OPT} -o${RAW_FIX_LASEPARATE_OLEN}"
    fi
    if [[ -n ${RAW_FIX_LASEPARATE_RLEN} && ${RAW_FIX_LASEPARATE_RLEN} -gt 0 ]]
    then 
        FIX_LASEPARATE_OPT="${FIX_LASEPARATE_OPT} -l${RAW_FIX_LASEPARATE_RLEN}"
    fi 
    if [[ -n ${RAW_FIX_LASEPARATE_USEREPEAT} ]]
    then 
    	ptype=""
    	if [[ "$1" -eq 0 ]]
    	then 
    		ptype="dalign"
    	elif [[ "$1" -eq 1 ]]
    	then 
    		ptype="repcomp"
    	fi
    	RAW_FIX_LASEPARATE_REPEAT="repeats_c${RAW_COV}_l${RAW_FIX_LAREPEAT_LEAVE_COV}h${RAW_FIX_LAREPEAT_ENTER_COV}_${ptype}_${RAW_REPMASK_LAREPEAT_REPEATTRACK}_${RAW_REPMASK_TANMASK_TRACK}_dust"
        FIX_LASEPARATE_OPT="${FIX_LASEPARATE_OPT} -r${RAW_FIX_LASEPARATE_REPEAT}"
    fi 

    # type is passed as argument
    FIX_LASEPARATE_OPT="${FIX_LASEPARATE_OPT} -T$1"
}

function setlassortOptions()
{
	FIX_LASSORT_OPT=""
	
	if [[ -z ${RAW_FIX_LASSORT_THREADS} ]]
	then 
		RAW_FIX_LASSORT_THREADS=8
	fi	
	FIX_LASSORT_OPT="${FIX_LASSORT_OPT} -t${RAW_FIX_LASSORT_THREADS}"
	
	if [[ -z ${RAW_FIX_LASSORT_MERGEFAN} ]]
	then 
		RAW_FIX_LASSORT_MERGEFAN=64
	fi	
	FIX_LASSORT_OPT="${FIX_LASSORT_OPT} -f${RAW_FIX_LASSORT_MERGEFAN}"

	if [[ -z ${RAW_FIX_LASSORT_SORT} ]]
	then 
		RAW_FIX_LASSORT_SORT=full
	fi	
	FIX_LASSORT_OPT="${FIX_LASSORT_OPT} -s${RAW_FIX_LASSORT_SORT}"
}

nblocks=$(getNumOfDbBlocks)

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

if [[ -z ${RAW_DALIGN_OUTDIR} ]]
then
	RAW_DALIGN_OUTDIR="dalign"
fi
if [[ -z ${RAW_REPCOMP_OUTDIR} ]]
then
	RAW_REPCOMP_OUTDIR="repcomp"
fi
if [[ -z ${RAW_DACCORD_OUTDIR} ]]
then
	RAW_DACCORD_OUTDIR="daccord"
fi

myTypes=("01-createSubdir, 02-daligner, 03-LAmerge, 04-LArepeat, 05-TKmerge, 06-TKcombine, 07-LAfilter, 08-LAq, 09-TKmerge, 10-LAfix" \
"01-createSubdir, 02-LAseparate, 03-repcomp, 04-LAmerge, 05-LArepeat, 06-TKmerge, 07-TKcombine, 08-LAq, 09-TKmerge, 10-LAfix" \
"01-createSubdir, 02-lassort2, 03-computeIntrinsicQV, 04_Catrack, 05_lasdetectsimplerepeats, 06_mergeAndSortRepeats, 07_lasfilteralignments, 08_mergesym2, 09_filtersym, 10_lasfilteralignmentsborderrepeats, 11_mergesym2, 12_filtersym, 13_filterchainsraw, 14_LAfilterChains, 15_LAfilter, 16_split, 16_LAmerge, 17_LAfix" \
"01_patchStats")
#type-0 - steps[1-10]: 01-createSubdir, 02-daligner, 03-LAmerge, 04-LArepeat, 05-TKmerge, 06-TKcombine, 07-LAfilter, 08-LAq, 09-TKmerge, 10-LAfix
#type-1 - steps[1-10]: 01-createSubdir, 02-LAseparate, 03-repcomp, 04-LAmerge, 05-LArepeat, 06-TKmerge, 07-TKcombine, 08-LAq, 09-TKmerge, 10-LAfix
#type-2 - steps[1-17]: 01-createSubdir, 02-lassort2, 03-computeIntrinsicQV, 04_Catrack, 05_lasdetectsimplerepeats, 06_mergeAndSortRepeats, 07_lasfilteralignments, 08_mergesym2, 09_filtersym, 10_lasfilteralignmentsborderrepeats, 11_mergesym2, 12_filtersym, 13_filterchainsraw, 14_LAfilterChains, 15_LAfilter, 16_split, 16_LAmerge, 17_LAfix
#type-3 - steps[1-1]:  01_patchStats 
if [[ ${RAW_PATCH_TYPE} -eq 0 ]]
then 
	if [[ ${currentStep} -eq 1 ]]
    then
        ### clean up plans 
        for x in $(ls fix_01_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        
        echo "if [[ -d ${RAW_DALIGN_OUTDIR} ]]; then mv ${RAW_DALIGN_OUTDIR} ${RAW_DALIGN_OUTDIR}_$(date '+%Y-%m-%d_%H-%M-%S'); fi && mkdir ${RAW_DALIGN_OUTDIR} && ln -s -r .${RAW_DB%.db}.* ${RAW_DB%.db}.db .${RAW_DAZZ_DB%.db}.* ${RAW_DAZZ_DB%.db}.db ${RAW_DALIGN_OUTDIR}" > fix_01_createSubdir_single_${RAW_DB%.db}.${slurmID}.plan
		for x in $(seq 1 ${nblocks})
	        do
			echo "mkdir -p ${RAW_DALIGN_OUTDIR}/d${x}"
		done >> fix_01_createSubdir_single_${RAW_DB%.db}.${slurmID}.plan
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_01_createSubdir_single_${RAW_DB%.db}.${slurmID}.version
    elif [[ ${currentStep} -eq 2 ]]
    then
        ### clean up plans 
        for x in $(ls fix_02_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        
        myCWD=$(pwd)
        ### find and set daligner options 
        setDalignerOptions
        ### create daligner commands
        cmdLine=1
        for x in $(seq 1 ${nblocks})
        do 

            if [[ -n ${RAW_FIX_DALIGNER_NUMACTL} && ${RAW_FIX_DALIGNER_NUMACTL} -gt 0 ]] && [[ "x${SLURM_NUMACTL}" == "x" || ${SLURM_NUMACTL} -eq 0 ]]
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
            cmd="cd ${RAW_DALIGN_OUTDIR} && PATH=${DAZZLER_PATH}/bin:\${PATH} ${NUMACTL}${DAZZLER_PATH}/bin/daligner${FIX_DALIGNER_OPT} ${RAW_DAZZ_DB%.db}.${x} ${RAW_DAZZ_DB%.db}.@${x}"
            cmdLine=$((${cmdLine}+1))
            count=0

            for y in $(seq ${x} ${nblocks})
            do  
                if [[ $count -lt ${RAW_FIX_DALIGNER_DAL} ]]
                then
                    count=$((${count}+1))
                else    
                    echo -n "${cmd}-$((y-1)) && mv"
                    z=${count}
		    while [[ $z -ge 1 ]]
		    do
			echo -n " ${RAW_DAZZ_DB%.db}.${x}.${RAW_DAZZ_DB%.db}.$((y-z)).las"
			z=$((z-1))
		    done
		    echo -n " d${x}"
		    if [[ -z "${RAW_FIX_DALIGNER_ASYMMETRIC}" ]]
		    then
			z=${count}
	                while [[ $z -ge 1 ]]
        	        do
				if [[ ${x} -ne $((y-z)) ]]
				then
                       		   echo -n " && mv ${RAW_DAZZ_DB%.db}.$((y-z)).${RAW_DAZZ_DB%.db}.${x}.las d$((y-z))"
				fi
                        	z=$((z-1)) 
                    	done   
		    fi
		    echo " && cd ${myCWD}"
                    if [[ -n ${RAW_FIX_DALIGNER_NUMACTL} && ${RAW_FIX_DALIGNER_NUMACTL} -gt 0 ]] && [[ "x${SLURM_NUMACTL}" == "x" || ${SLURM_NUMACTL} -eq 0 ]]
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
                    cmd="cd ${RAW_DALIGN_OUTDIR} && PATH=${DAZZLER_PATH}/bin:\${PATH} ${NUMACTL}${DAZZLER_PATH}/bin/daligner${FIX_DALIGNER_OPT} ${RAW_DAZZ_DB%.db}.${x} ${RAW_DAZZ_DB%.db}.@${y}"
                    cmdLine=$((${cmdLine}+1))
                    count=1
                fi
            done
            echo -n "${cmd}-${y} && mv"
            z=$((count-1))
                    while [[ $z -ge 0 ]]
                    do
                        echo -n " ${RAW_DAZZ_DB%.db}.${x}.${RAW_DAZZ_DB%.db}.$((y-z)).las"
                        z=$((z-1))
                    done
                    echo -n " d${x}"
                    if [[ -z "${RAW_FIX_DALIGNER_ASYMMETRIC}" ]]
                    then
                        z=$((count-1))
                        while [[ $z -ge 0 ]]
                        do
                                if [[ ${x} -ne $((y-z)) ]]
                                then
                                   echo -n " && mv ${RAW_DAZZ_DB%.db}.$((y-z)).${RAW_DAZZ_DB%.db}.${x}.las d$((y-z))"
                                fi
                                z=$((z-1))
                        done
                    fi
                    echo " && cd ${myCWD}"
    	done > fix_02_daligner_block_${RAW_DB%.db}.${slurmID}.plan
        echo "DAZZLER daligner $(git --git-dir=${DAZZLER_SOURCE_PATH}/DALIGNER/.git rev-parse --short HEAD)" > fix_02_daligner_block_${RAW_DB%.db}.${slurmID}.version
    elif [[ ${currentStep} -eq 3 ]]
    then
        ### clean up plans 
        for x in $(ls fix_03_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        myCWD=$(pwd)
        ### find and set LAmerge options 
        setLAmergeOptions
        ### create LAmerge commands
        for x in $(seq 1 ${nblocks})
        do 
            echo "cd ${RAW_DALIGN_OUTDIR} && ${MARVEL_PATH}/bin/LAmerge${FIX_LAMERGE_OPT} ${RAW_DB%.db} ${RAW_DAZZ_DB%.db}.dalign.${x}.las d${x} && cd ${myCWD}"
    	done > fix_03_LAmerge_block_${RAW_DB%.db}.${slurmID}.plan  
        echo "MARVEL LAmerge $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_03_LAmerge_block_${RAW_DB%.db}.${slurmID}.version       
    elif [[ ${currentStep} -eq 4 ]]
    then
        ### clean up plans 
        for x in $(ls fix_04_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        myCWD=$(pwd)
        setLArepeatOptions 1
        ### create LArepeat commands
        for x in $(seq 1 ${nblocks})
        do 
            echo "cd ${RAW_DALIGN_OUTDIR} && ${MARVEL_PATH}/bin/LArepeat${FIX_LAREPEAT_OPT} -b ${x} ${RAW_DB%.db} ${RAW_DAZZ_DB%.db}.dalign.${x}.las && cd ${myCWD}/"
            echo "cd ${RAW_DALIGN_OUTDIR} && ${DAZZLER_PATH}/bin/REPmask -v -c${RAW_DAZZ_FIX_LAREPEAT_THRESHOLD} -n${RAW_DAZZ_FIX_LAREPEAT_REPEATTRACK} ${RAW_DAZZ_DB%.db} ${RAW_DAZZ_DB%.db}.dalign.${x}.las && cd ${myCWD}/"
    	done > fix_04_LArepeat_block_${RAW_DB%.db}.${slurmID}.plan 
        echo "MARVEL LArepeat $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_04_LArepeat_block_${RAW_DB%.db}.${slurmID}.version
        echo "DAZZLER REPmask $(git --git-dir=${DAZZLER_SOURCE_PATH}/DAMASKER/.git rev-parse --short HEAD)" >> fix_04_LArepeat_block_${RAW_DB%.db}.${slurmID}.version        
    elif [[ ${currentStep} -eq 5 ]]
    then
        ### clean up plans 
        for x in $(ls fix_05_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        myCWD=$(pwd)
        # we need the name of the repeat track, especially if the plan starts with step4
        setLArepeatOptions 1
        ### find and set TKmerge options 
        if [[ -z "${FIX_TKMERGE_OPT}" ]]
        then 
            setTKmergeOptions
        fi
               
        ### create TKmerge command
        echo "cd ${RAW_DALIGN_OUTDIR} && ${MARVEL_PATH}/bin/TKmerge${FIX_TKMERGE_OPT} ${RAW_DB%.db} ${RAW_FIX_LAREPEAT_REPEATTRACK} && cp .${RAW_DB%.db}.${RAW_FIX_LAREPEAT_REPEATTRACK}.a2 .${RAW_DB%.db}.${RAW_FIX_LAREPEAT_REPEATTRACK}.d2 ${myCWD} && cd ${myCWD}" > fix_05_TKmerge_single_${RAW_DB%.db}.${slurmID}.plan      
        echo "cd ${RAW_DALIGN_OUTDIR} && ${DAZZLER_PATH}/bin/Catrack${FIX_TKMERGE_OPT} -f -v ${RAW_DAZZ_DB%.db} ${RAW_DAZZ_FIX_LAREPEAT_REPEATTRACK} && cp .${RAW_DAZZ_DB%.db}.${RAW_DAZZ_FIX_LAREPEAT_REPEATTRACK}.anno .${RAW_DAZZ_DB%.db}.${RAW_DAZZ_FIX_LAREPEAT_REPEATTRACK}.data ${myCWD}/ && cd ${myCWD}/" >> fix_05_TKmerge_single_${RAW_DB%.db}.${slurmID}.plan
        echo "MARVEL TKmerge $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_05_TKmerge_single_${RAW_DB%.db}.${slurmID}.version
        echo "DAZZLER Catrack $(git --git-dir=${DAZZLER_SOURCE_PATH}/DAZZ_DB/.git rev-parse --short HEAD)" >> fix_05_TKmerge_single_${RAW_DB%.db}.${slurmID}.version   
    elif [[ ${currentStep} -eq 6 ]]
    then
        ### clean up plans 
        for x in $(ls fix_06_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done     
        # we need the name of the repeat track, especially if the plan starts with step5
        setLArepeatOptions 1
        ### find and set TKcombine options
        setTKcombineOptions 1
        ### set repmask tracks 
        if [[ ${#RAW_REPMASK_LAREPEAT_COV[*]} -ne ${#RAW_REPMASK_BLOCKCMP[*]} ]]
        then 
            (>&2 echo "step ${currentStep} in RAW_PATCH_TYPE ${RAW_PATCH_TYPE}: arrays RAW_REPMASK_LAREPEAT_COV and RAW_REPMASK_BLOCKCMP must have same number of elements")
            exit 1
        fi
        RAW_REPMASK_REPEATTRACK=""
        for x in $(seq 1 ${#RAW_REPMASK_BLOCKCMP[*]})
        do
            idx=$(($x-1))
            RAW_REPMASK_REPEATTRACK="${RAW_REPMASK_REPEATTRACK} ${RAW_REPMASK_LAREPEAT_REPEATTRACK}_B${RAW_REPMASK_BLOCKCMP[${idx}]}C${RAW_REPMASK_LAREPEAT_COV[${idx}]}"
        done 
        ### create TKcombine command        
        if [[ -n ${RAW_REPMASK_REPEATTRACK} ]]
        then
            echo "${MARVEL_PATH}/bin/TKcombine${FIX_TKCOMBINE_OPT} ${RAW_DB%.db} ${RAW_FIX_LAREPEAT_REPEATTRACK}_${RAW_REPMASK_LAREPEAT_REPEATTRACK} ${RAW_FIX_LAREPEAT_REPEATTRACK} ${RAW_REPMASK_REPEATTRACK}" > fix_06_TKcombine_single_${RAW_DB%.db}.${slurmID}.plan
            echo "${MARVEL_PATH}/bin/TKcombine${FIX_TKCOMBINE_OPT} ${RAW_DB%.db} ${RAW_FIX_LAREPEAT_REPEATTRACK}_${RAW_REPMASK_LAREPEAT_REPEATTRACK}_${RAW_REPMASK_TANMASK_TRACK}_dust ${RAW_FIX_LAREPEAT_REPEATTRACK}_${RAW_REPMASK_LAREPEAT_REPEATTRACK} ${RAW_REPMASK_TANMASK_TRACK} dust" >> fix_06_TKcombine_single_${RAW_DB%.db}.${slurmID}.plan         
        else
            echo "ln -s .${RAW_DB%.db}.${RAW_FIX_LAREPEAT_REPEATTRACK}.d2 .${RAW_DB%.db}.${RAW_FIX_LAREPEAT_REPEATTRACK}_${RAW_REPMASK_LAREPEAT_REPEATTRACK}.d2"  > fix_06_TKcombine_single_${RAW_DB%.db}.${slurmID}.plan         
            echo "ln -s .${RAW_DB%.db}.${RAW_FIX_LAREPEAT_REPEATTRACK}.a2 .${RAW_DB%.db}.${RAW_FIX_LAREPEAT_REPEATTRACK}_${RAW_REPMASK_LAREPEAT_REPEATTRACK}.a2"  >> fix_06_TKcombine_single_${RAW_DB%.db}.${slurmID}.plan         
        fi 
        echo "MARVEL TKcombine $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_06_TKcombine_single_${RAW_DB%.db}.${slurmID}.version         
    elif [[ ${currentStep} -eq 7 ]]
    then
        ### clean up plans 
        for x in $(ls fix_07_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 

        mkdir -p identity
        ### create LAfilter commands - filter out identity overlaps - has to be done because revcomp and forcealign will loose those 
        for x in $(seq 1 ${nblocks})
        do  
            echo "${MARVEL_PATH}/bin/LAfilter -p -R 3 ${RAW_DB%.db} ${RAW_DALIGN_OUTDIR}/d${x}/${RAW_DAZZ_DB%.db}.${x}.${RAW_DAZZ_DB%.db}.${x}.las identity/${RAW_DAZZ_DB%.db}.identity.${x}.las"
		done > fix_07_LAfilter_block_${RAW_DB%.db}.${slurmID}.plan   
    echo "MARVEL LAfilter $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_07_LAfilter_block_${RAW_DB%.db}.${slurmID}.version    
    elif [[ ${currentStep} -eq 8 ]]
    then
        ### clean up plans 
        for x in $(ls fix_08_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        myCWD=$(pwd)
        ### find and set LAq options 
        setLAqOptions
        ### create LAq commands
        for x in $(seq 1 ${nblocks})
        do 
            echo "cd ${RAW_DALIGN_OUTDIR} && ${MARVEL_PATH}/bin/LAq${FIX_LAQ_OPT} -T trim0_d${RAW_FIX_LAQ_QTRIMCUTOFF}_s${RAW_FIX_LAQ_MINSEG}_dalign -Q q0_d${RAW_FIX_LAQ_QTRIMCUTOFF}_s${RAW_FIX_LAQ_MINSEG}_dalign ${RAW_DB%.db} -b ${x} ${RAW_DALIGN_OUTDIR}/${RAW_DAZZ_DB%.db}.dalign.${x}.las && cd ${myCWD}"
    	done > fix_08_LAq_block_${RAW_DB%.db}.${slurmID}.plan 
        echo "MARVEL LAq $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_08_LAq_block_${RAW_DB%.db}.${slurmID}.version                
    elif [[ ${currentStep} -eq 9 ]]
    then
        ### clean up plans 
        for x in $(ls fix_09_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        # we need the name of the q and trim track names, especially if the plan starts with step11
        if [[ -z ${FIX_LAREPEAT_OPT} ]]
        then 
            setLAqOptions
        fi  
        if [[ -z ${FIX_TKMERGE_OPT} ]]
        then 
            setTKmergeOptions
        fi
        myCWD=$(pwd)
        ### create TKmerge command
        echo "cd ${RAW_DALIGN_OUTDIR} && ${MARVEL_PATH}/bin/TKmerge${FIX_TKMERGE_OPT} ${RAW_DB%.db} trim0_d${RAW_FIX_LAQ_QTRIMCUTOFF}_s${RAW_FIX_LAQ_MINSEG}_dalign && cp .${RAW_DB%.db}.trim0_d${RAW_FIX_LAQ_QTRIMCUTOFF}_s${RAW_FIX_LAQ_MINSEG}_dalign.a2 .${RAW_DB%.db}.trim0_d${RAW_FIX_LAQ_QTRIMCUTOFF}_s${RAW_FIX_LAQ_MINSEG}_dalign.d2 ${myCWD}/ && cd ${myCWD}" > fix_09_TKmerge_block_${RAW_DB%.db}.${slurmID}.plan
        echo "cd ${RAW_DALIGN_OUTDIR} && ${MARVEL_PATH}/bin/TKmerge${FIX_TKMERGE_OPT} ${RAW_DB%.db} q0_d${RAW_FIX_LAQ_QTRIMCUTOFF}_s${RAW_FIX_LAQ_MINSEG}_dalign && cp .${RAW_DB%.db}.q0_d${RAW_FIX_LAQ_QTRIMCUTOFF}_s${RAW_FIX_LAQ_MINSEG}_dalign.a2 .${RAW_DB%.db}.q0_d${RAW_FIX_LAQ_QTRIMCUTOFF}_s${RAW_FIX_LAQ_MINSEG}_dalign.d2 ${myCWD}/ && cd ${myCWD}" >> fix_09_TKmerge_block_${RAW_DB%.db}.${slurmID}.plan       
        echo "MARVEL TKmerge $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_09_TKmerge_block_${RAW_DB%.db}.${slurmID}.version               
	### 10_LAfix    
    elif [[ ${currentStep} -eq 10 ]]
    then
        ### clean up plans 
        for x in $(ls fix_10_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        ### find and set LAfix options 
        setLAfixOptions dalign
        mkdir -p ${RAW_FIX_LAFIX_PATH}
		
		addopt=""

        for x in $(seq 1 ${nblocks})
        do 
        	if [[ -n ${RAW_FIX_LAFIX_TRIMFILEPREFIX} ]]
        	then 
        		addopt="-T${RAW_FIX_LAFIX_TRIMFILEPREFIX}_${x}.txt "
        	fi
            echo "${MARVEL_PATH}/bin/LAfix${FIX_LAFIX_OPT} ${addopt}${RAW_DB%.db} ${RAW_DALIGN_OUTDIR}/${RAW_DAZZ_DB%.db}.dalign.${x}.las ${RAW_FIX_LAFIX_PATH}/${RAW_DB%.db}.${x}${RAW_FIX_LAFIX_FILESUFFIX}.fasta"
    	done > fix_10_LAfix_block_${RAW_DB%.db}.${slurmID}.plan
    echo "MARVEL LAfix $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_10_LAfix_block_${RAW_DB%.db}.${slurmID}.version                 
    else 
        (>&2 echo "step ${currentStep} in RAW_PATCH_TYPE ${RAW_PATCH_TYPE} not supported")
        (>&2 echo "valid steps are: ${myTypes[${RAW_PATCH_TYPE}]}")
        exit 1        
    fi     
#type-1 steps   [ 1-1] :  01-createSubdir, 02-LAseparate, 03-repcomp, 04-LAmerge, 05-LArepeat, 06-TKmerge, 07-TKcombine, 08-LAq, 09-TKmerge, 10-LAfix     
elif [[ ${RAW_PATCH_TYPE} -eq 1 ]]
then
	if [[ ${currentStep} -eq 1 ]]
    then
		### clean up plans 
	    for x in $(ls fix_01_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
	    do            
	        rm $x
	    done 
	    
	    echo "if [[ -d ${RAW_REPCOMP_OUTDIR} ]]; then mv ${RAW_REPCOMP_OUTDIR} ${RAW_REPCOMP_OUTDIR}_$(date '+%Y-%m-%d_%H-%M-%S'); fi && mkdir ${RAW_REPCOMP_OUTDIR} && ln -s -r .${RAW_DB%.db}.* ${RAW_DB%.db}.db .${RAW_DAZZ_DB%.db}.* ${RAW_DAZZ_DB%.db}.db ${RAW_REPCOMP_OUTDIR}" > fix_01_createSubdir_single_${RAW_DB%.db}.${slurmID}.plan
		for x in $(seq 1 ${nblocks})
	    do
			echo "mkdir -p ${RAW_REPCOMP_OUTDIR}/r${x}"
			echo "mkdir -p ${RAW_REPCOMP_OUTDIR}/d${x}_ForRepComp"
			echo "mkdir -p ${RAW_REPCOMP_OUTDIR}/d${x}_NoRepComp"
		done >> fix_01_createSubdir_single_${RAW_DB%.db}.${slurmID}.plan
	    echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_01_createSubdir_single_${RAW_DB%.db}.${slurmID}.version
  	#### LAseparate
  	elif [[ ${currentStep} -eq 2 ]]
    then
        ### clean up plans 
        for x in $(ls fix_02_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 

        ### find and set LAseparate options 
        setLAseparateOptions 0

        for x in $(seq 1 ${nblocks}); 
        do 
            for y in $(seq 1 ${nblocks}); 
            do 
                if [[ ! -f ${RAW_DALIGN_OUTDIR}/d${x}/${RAW_DAZZ_DB%.db}.${x}.${RAW_DAZZ_DB%.db}.${y}.las ]]
                then
                    (>&2 echo "step ${currentStep} in RAW_PATCH_TYPE ${RAW_PATCH_TYPE}: File missing ${RAW_DALIGN_OUTDIR}/d${x}/${RAW_DAZZ_DB%.db}.${x}.${RAW_DAZZ_DB%.db}.${y}.las!!")
                    exit 1                    
                fi
                echo "${MARVEL_PATH}/bin/LAseparate${FIX_LASEPARATE_OPT} ${RAW_DB%.db} ${RAW_DALIGN_OUTDIR}/d${x}/${RAW_DAZZ_DB%.db}.${x}.${RAW_DAZZ_DB%.db}.${y}.las ${RAW_REPCOMP_OUTDIR}/d${x}_ForRepComp/${RAW_DAZZ_DB%.db}.${x}.${RAW_DAZZ_DB%.db}.${y}.las ${RAW_REPCOMP_OUTDIR}/d${x}_ForRepComp/${RAW_DAZZ_DB%.db}.${x}.${RAW_DAZZ_DB%.db}.${y}.las"                
            done 
    	done > fix_02_LAseparate_block_${RAW_DB%.db}.${slurmID}.plan
    	echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_02_LAseparate_block_${RAW_DB%.db}.${slurmID}.version
    #### repcomp 
    elif [[ ${currentStep} -eq 3 ]]
    then
        ### clean up plans 
        for x in $(ls fix_03_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        ### find and set repcomp options 
        setRepcompOptions

        cmdLine=1
        for x in $(seq 1 ${nblocks}); 
        do 
            srcDir=${RAW_REPCOMP_OUTDIR}/d${x}_ForRepComp
            desDir=${RAW_REPCOMP_OUTDIR}/r${x}

            if [[ ! -d ${desDir} ]]
            then
                mkdir -p ${desDir}
            fi
            start=${x}

            for y in $(seq ${start} ${nblocks}); 
            do 
                movDir=${RAW_REPCOMP_OUTDIR}/r${y}
                if [[ -f ${srcDir}/${RAW_DAZZ_DB%.db}.${x}.${RAW_DAZZ_DB%.db}.${y}.las ]]
                then 
                    if [[ -n ${RAW_FIX_REPCOMP_NUMACTL} && ${RAW_FIX_REPCOMP_NUMACTL} -gt 0 ]] && [[ "x${SLURM_NUMACTL}" == "x" || ${SLURM_NUMACTL} -eq 0 ]]
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
                    echo -n "${NUMACTL}${REPCOMP_PATH}/bin/repcomp${FIX_REPCOMP_OPT} -T/tmp/${RAW_DAZZ_DB%.db}.${x}.${y} ${desDir}/${RAW_DAZZ_DB%.db}.repcomp.${x}.${y} ${RAW_DAZZ_DB%.db} ${srcDir}/${RAW_DAZZ_DB%.db}.${x}.${RAW_DAZZ_DB%.db}.${y}.las"
                    cmdLine=$((${cmdLine}+1))
                    if [[ $x -eq $y ]]
                    then
                        echo ""
                    else    
                        echo " && mv ${desDir}/${RAW_DAZZ_DB%.db}.repcomp.${x}.${y}_r.las ${movDir}/"
                    fi
                else
                    (>&2 echo "step ${currentStep} in RAW_FIX_TYPE ${RAW_FIX_TYPE}: File missing ${srcDir}/${RAW_DAZZ_DB%.db}.${x}.${RAW_DAZZ_DB%.db}.${y}.las!!")
                    exit 1
                fi
            done 
		done > fix_03_repcomp_block_${RAW_DB%.db}.${slurmID}.plan
    	echo "repcomp $(git --git-dir=${REPCOMP_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_03_repcomp_block_${RAW_DB%.db}.${slurmID}.version
	### 04_LAmergeLAfilter
    elif [[ ${currentStep} -eq 4 ]]
    then
        ### clean up plans 
        for x in $(ls fix_04_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        myCWD=$(pwd)
        ### find and set LAmerge options 
        setLAmergeOptions
        setRepcompOptions
        ### create LAmerge commands
        for x in $(seq 1 ${nblocks})
        do 
            echo "cd ${RAW_REPCOMP_OUTDIR} && ${MARVEL_PATH}/bin/LAmerge${FIX_LAMERGE_OPT} ${RAW_DB%.db} ${RAW_DAZZ_DB%.db}.repcomp.${x}.las r${x} ${RAW_REPCOMP_OUTDIR}/d${x}_ForRepComp ${RAW_REPCOMP_OUTDIR}/d${x}_NoRepComp identity/${RAW_DAZZ_DB%.db}.identity.${x}.las && ${MARVEL_PATH}/bin/LAfilter -p -R6 ${RAW_DB%.db} ${RAW_DAZZ_DB%.db}.repcomp.${x}.las ${RAW_DAZZ_DB%.db}.repcompFilt.${x}.las && cd ${myCWD}"                                                                                                                     
    	done > fix_04_LAmerge_block_${RAW_DB%.db}.${slurmID}.plan
    	echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_04_LAmerge_block_${RAW_DB%.db}.${slurmID}.version      
    ### 05_LArepeat
    elif [[ ${currentStep} -eq 5 ]]
    then
        ### clean up plans 
        for x in $(ls fix_05_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        myCWD=$(pwd)
        setLArepeatOptions 2
        ### create LArepeat commands
        for x in $(seq 1 ${nblocks})
        do 
            echo "cd ${RAW_REPCOMP_OUTDIR} && ${MARVEL_PATH}/bin/LArepeat${FIX_LAREPEAT_OPT} -b ${x} ${RAW_DB%.db} ${RAW_DAZZ_DB%.db}.repcompFilt.${x}.las && cd ${myCWD}/"
            echo "cd ${RAW_REPCOMP_OUTDIR} && ${DAZZLER_PATH}/bin/REPmask -v -c${RAW_DAZZ_FIX_LAREPEAT_THRESHOLD} -n${RAW_DAZZ_FIX_LAREPEAT_REPEATTRACK} ${RAW_DAZZ_DB%.db} ${RAW_DAZZ_DB%.db}.repcompFilt.${x}.las && cd ${myCWD}/"
            
    	done > fix_05_LArepeat_block_${RAW_DB%.db}.${slurmID}.plan
    	echo "MARVEL LArepeat $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_05_LArepeat_block_${RAW_DB%.db}.${slurmID}.version
    	echo "DAZZLER REPmask $(git --git-dir=${DAZZLER_SOURCE_PATH}/DAMASKER/.git rev-parse --short HEAD)" >> fix_05_LArepeat_block_${RAW_DB%.db}.${slurmID}.version
    ### 06_TKmerge         
    elif [[ ${currentStep} -eq 6 ]]
    then
        ### clean up plans 
        for x in $(ls fix_06_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        myCWD=$(pwd)
        # we need the name of the repeat track, especially if the plan starts with step4
        setLArepeatOptions 2
        ### find and set TKmerge options 
        if [[ -z ${FIX_TKMERGE_OPT} ]]
        then 
            setTKmergeOptions
        fi
        ### create TKmerge command
        echo "cd ${RAW_REPCOMP_OUTDIR} && ${MARVEL_PATH}/bin/TKmerge${FIX_TKMERGE_OPT} ${RAW_DB%.db} ${RAW_FIX_LAREPEAT_REPEATTRACK} && cp .${RAW_DB%.db}.${RAW_FIX_LAREPEAT_REPEATTRACK}.a2 .${RAW_DB%.db}.${RAW_FIX_LAREPEAT_REPEATTRACK}.d2 ${myCWD} && cd ${myCWD}" > fix_06_TKmerge_single_${RAW_DB%.db}.${slurmID}.plan      
        echo "cd ${RAW_REPCOMP_OUTDIR} && ${DAZZLER_PATH}/bin/Catrack${FIX_TKMERGE_OPT} -f -v ${RAW_DAZZ_DB%.db} ${RAW_DAZZ_FIX_LAREPEAT_REPEATTRACK} && cp .${RAW_DAZZ_DB%.db}.${RAW_DAZZ_FIX_LAREPEAT_REPEATTRACK}.anno .${RAW_DAZZ_DB%.db}.${RAW_DAZZ_FIX_LAREPEAT_REPEATTRACK}.data ${myCWD}/ && cd ${myCWD}/" >> fix_06_TKmerge_single_${RAW_DB%.db}.${slurmID}.plan
        echo "MARVEL TKmerge $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_06_TKmerge_single_${RAW_DB%.db}.${slurmID}.version
        echo "DAZZLER Catrack $(git --git-dir=${DAZZLER_SOURCE_PATH}/DAZZ_DB/.git rev-parse --short HEAD)" >> fix_06_TKmerge_single_${RAW_DB%.db}.${slurmID}.version
    ### 07_TKcombine   
    elif [[ ${currentStep} -eq 7 ]]
    then
        ### clean up plans 
        for x in $(ls fix_07_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done     
        # we need the name of the repeat track, especially if the plan starts with step5
        setLArepeatOptions 2
        ### find and set TKcombine options
        setTKcombineOptions 1
        ### set repmask tracks 
        if [[ ${#RAW_REPMASK_LAREPEAT_COV[*]} -ne ${#RAW_REPMASK_BLOCKCMP[*]} ]]
        then 
            (>&2 echo "step ${currentStep} in RAW_PATCH_TYPE ${RAW_PATCH_TYPE}: arrays RAW_REPMASK_LAREPEAT_COV and RAW_REPMASK_BLOCKCMP must have same number of elements")
            exit 1
        fi
        RAW_REPMASK_REPEATTRACK=""
        for x in $(seq 1 ${#RAW_REPMASK_BLOCKCMP[*]})
        do
            idx=$(($x-1))
            RAW_REPMASK_REPEATTRACK="${RAW_REPMASK_REPEATTRACK} ${RAW_REPMASK_LAREPEAT_REPEATTRACK}_B${RAW_REPMASK_BLOCKCMP[${idx}]}C${RAW_REPMASK_LAREPEAT_COV[${idx}]}"
        done 
        ### create TKcombine command        
        if [[ -n ${RAW_REPMASK_REPEATTRACK} ]]
        then
            echo "${MARVEL_PATH}/bin/TKcombine${FIX_TKCOMBINE_OPT} ${RAW_DB%.db} ${RAW_FIX_LAREPEAT_REPEATTRACK}_${RAW_REPMASK_LAREPEAT_REPEATTRACK} ${RAW_FIX_LAREPEAT_REPEATTRACK} ${RAW_REPMASK_REPEATTRACK}" > fix_07_TKcombine_single_${RAW_DB%.db}.${slurmID}.plan
            echo "${MARVEL_PATH}/bin/TKcombine${FIX_TKCOMBINE_OPT} ${RAW_DB%.db} ${RAW_FIX_LAREPEAT_REPEATTRACK}_${RAW_REPMASK_LAREPEAT_REPEATTRACK}_${RAW_REPMASK_TANMASK_TRACK}_dust ${RAW_FIX_LAREPEAT_REPEATTRACK}_${RAW_REPMASK_LAREPEAT_REPEATTRACK} ${RAW_REPMASK_TANMASK_TRACK} dust" >> fix_07_TKcombine_single_${RAW_DB%.db}.${slurmID}.plan         
        else
            echo "ln -s .${RAW_DB%.db}.${RAW_FIX_LAREPEAT_REPEATTRACK}.d2 .${RAW_DB%.db}.${RAW_FIX_LAREPEAT_REPEATTRACK}_${RAW_REPMASK_LAREPEAT_REPEATTRACK}.d2"  > fix_07_TKcombine_single_${RAW_DB%.db}.${slurmID}.plan         
            echo "ln -s .${RAW_DB%.db}.${RAW_FIX_LAREPEAT_REPEATTRACK}.a2 .${RAW_DB%.db}.${RAW_FIX_LAREPEAT_REPEATTRACK}_${RAW_REPMASK_LAREPEAT_REPEATTRACK}.a2"  >> fix_07_TKcombine_single_${RAW_DB%.db}.${slurmID}.plan         
        fi 
        echo "MARVEL TKcombine $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_07_TKcombine_single_${RAW_DB%.db}.${slurmID}.version
    ### 08_LAq  
    elif [[ ${currentStep} -eq 8 ]]
    then
        ### clean up plans 
        for x in $(ls fix_08_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        myCWD=$(pwd)
        ### find and set LAq options 
        setLAqOptions
        ### create LAq commands
        for x in $(seq 1 ${nblocks})
        do 
            echo "cd ${RAW_REPCOMP_OUTDIR} && ${MARVEL_PATH}/bin/LAq${FIX_LAQ_OPT} -T trim0_d${RAW_FIX_LAQ_QTRIMCUTOFF}_s${RAW_FIX_LAQ_MINSEG}_repcomp -Q q0_d${RAW_FIX_LAQ_QTRIMCUTOFF}_s${RAW_FIX_LAQ_MINSEG}_repcomp ${RAW_DB%.db} -b ${x} ${RAW_DB%.db}.repcompFilt.${x}.las && cd ${myCWD}"
		done > fix_08_LAq_block_${RAW_DB%.db}.${slurmID}.plan
    	echo "MARVEL LAq $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_08_LAq_block_${RAW_DB%.db}.${slurmID}.version
	### 09_TKmerge    	                 
    elif [[ ${currentStep} -eq 9 ]]
    then
        ### clean up plans 
        for x in $(ls fix_09_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done  
        myCWD=$(pwd)   
        # we need the name of the q and trim track names, especially if the plan starts with step11
        if [[ -z ${FIX_LAQ_OPT} ]]
        then 
            setLAqOptions
        fi  
        if [[ -z ${FIX_TKMERGE_OPT} ]]
        then 
            setTKmergeOptions
        fi
        ### create TKmerge command
        echo "cd ${RAW_REPCOMP_OUTDIR} && ${MARVEL_PATH}/bin/TKmerge${FIX_TKMERGE_OPT} ${RAW_DB%.db} trim0_d${RAW_FIX_LAQ_QTRIMCUTOFF}_s${RAW_FIX_LAQ_MINSEG}_repcomp && cp .${RAW_DB%.db}.trim0_d${RAW_FIX_LAQ_QTRIMCUTOFF}_s${RAW_FIX_LAQ_MINSEG}_repcomp.a2 .${RAW_DB%.db}.trim0_d${RAW_FIX_LAQ_QTRIMCUTOFF}_s${RAW_FIX_LAQ_MINSEG}_repcomp.d2 ${myCWD}/ && cd ${myCWD}" > fix_09_TKmerge_block_${RAW_DB%.db}.${slurmID}.plan
        echo "cd ${RAW_REPCOMP_OUTDIR} && ${MARVEL_PATH}/bin/TKmerge${FIX_TKMERGE_OPT} ${RAW_DB%.db} q0_d${RAW_FIX_LAQ_QTRIMCUTOFF}_s${RAW_FIX_LAQ_MINSEG}_repcomp&& cp .${RAW_DB%.db}.q0_d${RAW_FIX_LAQ_QTRIMCUTOFF}_s${RAW_FIX_LAQ_MINSEG}_repcomp.a2 .${RAW_DB%.db}.q0_d${RAW_FIX_LAQ_QTRIMCUTOFF}_s${RAW_FIX_LAQ_MINSEG}_repcomp.d2 ${myCWD}/ && cd ${myCWD}" >> fix_09_TKmerge_block_${RAW_DB%.db}.${slurmID}.plan
        echo "MARVEL TKmerge $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_09_TKmerge_block_${RAW_DB%.db}.${slurmID}.version               
   ### 10 LAfix
    elif [[ ${currentStep} -eq 10 ]]
    then
        ### clean up plans 
        for x in $(ls fix_10_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        ### find and set LAfix options 
        setLAfixOptions repcomp
        mkdir -p ${RAW_FIX_LAFIX_PATH}
		addopt=""
        for x in $(seq 1 ${nblocks})
        do 
        	if [[ -n ${RAW_FIX_LAFIX_TRIMFILEPREFIX} ]]
        	then 
        		addopt="-T${RAW_FIX_LAFIX_TRIMFILEPREFIX}_${x}.txt "
        	fi
            echo "${MARVEL_PATH}/bin/LAfix${FIX_LAFIX_OPT} ${addopt}${RAW_DB%.db} ${RAW_REPCOMP_OUTDIR}/${RAW_DAZZ_DB%.db}.repcompFilt.${x}.las ${RAW_FIX_LAFIX_PATH}/${RAW_DB%.db}.${x}${RAW_FIX_LAFIX_FILESUFFIX}.fasta"
		done > fix_10_LAfix_block_${RAW_DB%.db}.${slurmID}.plan
    	echo "MARVEL LAfix $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_10_LAfix_block_${RAW_DB%.db}.${slurmID}.version                                  
    else 
        (>&2 echo "step ${currentStep} in RAW_PATCH_TYPE ${RAW_PATCH_TYPE} not supported")
        (>&2 echo "valid steps are: ${myTypes[${RAW_PATCH_TYPE}]}")
        exit 1        
    fi  
elif [[ ${RAW_PATCH_TYPE} -eq 3 ]]
then
  	if [[ ${currentStep} -eq 1 ]]
    then
        ### clean up plans 
        for x in $(ls fix_01_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done    
        
        if [[ -n ${SLURM_STATS} && ${SLURM_STATS} -gt 0 ]]
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
		if [[ -n ${MARVEL_STATS} && ${MARVEL_STATS} -gt 0 ]]
   		then
	        ### create assemblyStats plan
	    	echo "${SUBMIT_SCRIPTS_PATH}/patchingStats.sh ${configFile} 1" > fix_01_patchingStats_block_${FIX_DB%.db}.${slurmID}.plan
		fi
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_01_patchingStats_block_${FIX_DB%.db}.${slurmID}.version
    else
		(>&2 echo "step ${currentStep} in RAW_PATCH_TYPE ${RAW_PATCH_TYPE} not supported")
        (>&2 echo "valid steps are: ${myTypes[${RAW_PATCH_TYPE}]}")
        exit 1        
    fi
elif [[ ${RAW_PATCH_TYPE} -eq 2 ]]
then 
    ### create sub-directory and link relevant DB and Track files
    if [[ ${currentStep} -eq 1 ]]
    then
        ### clean up plans 
    for x in $(ls fix_01_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done                 

        echo "if [[ -d ${RAW_DACCORD_OUTDIR} ]]; then mv ${RAW_DACCORD_OUTDIR} ${RAW_DACCORD_OUTDIR}_$(date '+%Y-%m-%d_%H-%M-%S'); fi && mkdir ${RAW_DACCORD_OUTDIR} && ln -s -r .${RAW_DAZZ_DB%.db}.* ${RAW_DAZZ_DB%.db}.db ${RAW_DACCORD_OUTDIR}" > fix_01_patchingStats_single_${RAW_DB%.db}.${slurmID}.plan
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_01_patchingStats_single_${RAW_DB%.db}.${slurmID}.version
 	### 02-lassort
	elif [[ ${currentStep} -eq 2 ]]
    then
        ### clean up plans 
    	for x in $(ls fix_02_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 

		setlassortOptions
		
		for x in $(seq 1 ${nblocks})
        do
        	echo "${LASTOOLS_PATH}/bin/lassort ${FIX_LASSORT_OPT} ${RAW_DACCORD_OUTDIR}/${RAW_DB%.db}.${x}.${FIX_FILT_ENDING}sort.las ${FIX_DB%.db}.${x}.${FIX_FILT_ENDING}.las"
		done > filt_02_lassort_block_${FIX_DB%.db}.${slurmID}.plan    	         
        echo "LASTOOLS lassort $(git --git-dir=${DACCORD_SOURCE_PATH}/.git rev-parse --short HEAD)" > filt_02_lassort_block_${FIX_DB%.db}.${slurmID}.version
    ### 03-computeIntrinsicQV
	elif [[ ${currentStep} -eq 3 ]]
    then
        ### clean up plans 
        for x in $(ls filt_03_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 

		setLAfilterOptions
		
		for x in $(seq 1 ${fixblocks})
        do
        	echo "LIBMAUS2_DAZZLER_ALIGN_ALIGNMENTFILECONSTANTS_TRACE_XOVR=75 ${DACCORD_PATH}/bin/computeintrinsicqv2 -d${FIX_COV} ${FIX_FILT_OUTDIR}/${FIX_DAZZ_DB%.db}.db ${FIX_FILT_OUTDIR}/${FIX_DB%.db}.${x}.${FIX_FILT_ENDING}sort.las"
		done > filt_03_computeintrinsicqv2_block_${FIX_DB%.db}.${slurmID}.plan    	         
        echo "DACCORD computeintrinsicqv2 $(git --git-dir=${DACCORD_SOURCE_PATH}/.git rev-parse --short HEAD)" > filt_03_computeintrinsicqv2_block_${FIX_DB%.db}.${slurmID}.version
	### 04_Catrack
	elif [[ ${currentStep} -eq 4 ]]
    then
        ### clean up plans 
        for x in $(ls filt_04_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        echo "PATH=${DAZZLER_PATH}/bin:\${PATH} ${DAZZLER_PATH}/bin/Catrack -v -f -d ${FIX_FILT_OUTDIR}/${FIX_DAZZ_DB%.db}.db inqual" > filt_04_Catrack_single_${FIX_DB%.db}.${slurmID}.plan
		echo "DAZZ_DB Catrack $(git --git-dir=${DAZZLER_SOURCE_PATH}/DAZZ_DB/.git rev-parse --short HEAD)" > filt_04_Catrack_single_${FIX_DB%.db}.${slurmID}.version
                 
    ### 05_lasdetectsimplerepeats
    elif [[ ${currentStep} -eq 5 ]]
    then
        ### clean up plans 
        for x in $(ls filt_05_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        setLAfilterOptions
        
        OPT=""
        if [[ -z "${FIX_FILT_LASDETECTSIMPLEREPEATS_ERATE}" ]]
        then 
        	FIX_FILT_LASDETECTSIMPLEREPEATS_ERATE=0.35
   	 	fi 
   	 	
   	 	OPT="${OPT} -d$((FIX_COV/2)) -e${FIX_FILT_LASDETECTSIMPLEREPEATS_ERATE}"
    
        for x in $(seq 1 ${fixblocks})
        do
        	echo "LIBMAUS2_DAZZLER_ALIGN_ALIGNMENTFILECONSTANTS_TRACE_XOVR=75 ${DACCORD_PATH}/bin/lasdetectsimplerepeats ${OPT} ${FIX_FILT_OUTDIR}/${FIX_DAZZ_DB%.db}.${x}.rep ${FIX_FILT_OUTDIR}/${FIX_DAZZ_DB%.db}.db ${FIX_FILT_OUTDIR}/${FIX_DB%.db}.${x}.${FIX_FILT_ENDING}sort.las"
		done > filt_05_lasdetectsimplerepeats_block_${FIX_DB%.db}.${slurmID}.plan
      	echo "DACCORD lasdetectsimplerepeats $(git --git-dir=${DACCORD_SOURCE_PATH}/.git rev-parse --short HEAD)" > filt_05_lasdetectsimplerepeats_block_${FIX_DB%.db}.${slurmID}.version
    ### 06_mergeAndSortRepeats
    elif [[ ${currentStep} -eq 6 ]]
    then
        ### clean up plans 
        for x in $(ls filt_06_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
    	files="${FIX_FILT_OUTDIR}/${FIX_DAZZ_DB%.db}.[0-9].rep"
        if [[ ${fixblocks} -gt 9 ]]
        then
        	files="${files} ${FIX_FILT_OUTDIR}/${FIX_DAZZ_DB%.db}.[0-9][0-9].rep"
        elif [[ ${fixblocks} -gt 99 ]]
        then
        	files="${files} ${FIX_FILT_OUTDIR}/${FIX_DAZZ_DB%.db}.[0-9][0-9][0-9].rep"
        elif [[ ${fixblocks} -gt 999 ]]
        then
        	files="${files} ${FIX_FILT_OUTDIR}/${FIX_DAZZ_DB%.db}.[0-9][0-9][0-9][0-9].rep"
        elif [[ ${fixblocks} -gt 9999 ]]
        then
        	files="${files} ${FIX_FILT_OUTDIR}/${FIX_DAZZ_DB%.db}.[0-9][0-9][0-9][0-9][0-9].rep"
    	else
    		(>&2 echo "05_mergeAndSortRepeats: more than 99999 db blocks are not supported!!!")
        	exit 1	
    	fi
    
    	echo "cat ${files} > ${FIX_FILT_OUTDIR}/${FIX_DAZZ_DB%.db}.rep" > filt_06_mergeAndSortRepeats_single_${FIX_DB%.db}.${slurmID}.plan
    	echo "cat ${FIX_FILT_OUTDIR}/${FIX_DAZZ_DB%.db}.rep | ${DACCORD_PATH}/bin/repsort ${FIX_FILT_OUTDIR}/${FIX_DAZZ_DB%.db}.db > ${FIX_FILT_OUTDIR}/${FIX_DAZZ_DB%.db}.sort.rep" >> filt_06_mergeAndSortRepeats_single_${FIX_DB%.db}.${slurmID}.plan 
    	echo "rm ${files} ${FIX_FILT_OUTDIR}/${FIX_DAZZ_DB%.db}.rep" >> filt_06_mergeAndSortRepeats_single_${FIX_DB%.db}.${slurmID}.plan
        echo "DACCORD repsort $(git --git-dir=${DACCORD_SOURCE_PATH}/.git rev-parse --short HEAD)" > filt_06_mergeAndSortRepeats_single_${FIX_DB%.db}.${slurmID}.version
    ### 07_lasfilteralignments 
    elif [[ ${currentStep} -eq 7 ]]
    then
        ### clean up plans 
        for x in $(ls filt_07_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        setLAfilterOptions
        
        OPT=""
        
        if [[ -z "${FIX_FILT_LASFILTERALIGNMENTS_ERATE}" ]]
        then 
        	FIX_FILT_LASFILTERALIGNMENTS_ERATE=0.35
   	 	fi 
   	 	
   	 	OPT="${OPT} -e${FIX_FILT_LASFILTERALIGNMENTS_ERATE}"
    
        for x in $(seq 1 ${fixblocks})
        do
        	echo "LIBMAUS2_DAZZLER_ALIGN_ALIGNMENTFILECONSTANTS_TRACE_XOVR=75 ${DACCORD_PATH}/bin/lasfilteralignments ${OPT} ${FIX_FILT_OUTDIR}/${FIX_DB%.db}.${x}.${FIX_FILT_ENDING}LasFiltAln.las ${FIX_FILT_OUTDIR}/${FIX_DAZZ_DB%.db}.db ${FIX_FILT_OUTDIR}/${FIX_DB%.db}.${x}.${FIX_FILT_ENDING}sort.las"
		done > filt_07_lasfilteralignments_block_${FIX_DB%.db}.${slurmID}.plan
      	echo "DACCORD lasfilteralignments $(git --git-dir=${DACCORD_SOURCE_PATH}/.git rev-parse --short HEAD)" > filt_07_lasfilteralignments_block_${FIX_DB%.db}.${slurmID}.version
    ### 08_mergesym2
    elif [[ ${currentStep} -eq 8 ]]
    then
        ### clean up plans 
        for x in $(ls filt_08_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        setLAfilterOptions
        OPT=""
        echo "LIBMAUS2_DAZZLER_ALIGN_ALIGNMENTFILECONSTANTS_TRACE_XOVR=75 ${DACCORD_PATH}/bin/mergesym2 ${OPT} ${FIX_FILT_OUTDIR}/${FIX_DB%.db}.${FIX_FILT_ENDING}LasFiltAln.las.sym ${FIX_FILT_OUTDIR}/${FIX_DAZZ_DB%.db}.db ${FIX_FILT_OUTDIR}/${FIX_DB%.db}.*.${FIX_FILT_ENDING}LasFiltAln.las.sym" > filt_08_mergesym2_single_${FIX_DB%.db}.${slurmID}.plan
        echo "rm ${FIX_FILT_OUTDIR}/${FIX_DB%.db}.*.${FIX_FILT_ENDING}LasFiltAln.las.sym" >> filt_08_mergesym2_single_${FIX_DB%.db}.${slurmID}.plan
        echo "DACCORD mergesym2 $(git --git-dir=${DACCORD_SOURCE_PATH}/.git rev-parse --short HEAD)" > filt_08_mergesym2_single_${FIX_DB%.db}.${slurmID}.version        
	### 09_filtersym
    elif [[ ${currentStep} -eq 9 ]]
    then
        ### clean up plans 
        for x in $(ls filt_09_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        setLAfilterOptions
        OPT=""        
        
		if [[ -z "${FIX_FILT_FILTERSYM_VERBOSE}" ]]
        then
        	FIX_FILT_FILTERSYM_VERBOSE=1
   	 	fi 
   	 	
   	 	if [[ -n "${FIX_FILT_FILTERSYM_VERBOSE}" && ${FIX_FILT_FILTERSYM_VERBOSE} != 0 ]]
        then
   	 		OPT="--verbose" 
   	 	fi
   	 	
   	 	for x in $(seq 1 ${fixblocks})
        do
    		echo "LIBMAUS2_DAZZLER_ALIGN_ALIGNMENTFILECONSTANTS_TRACE_XOVR=75 ${DACCORD_PATH}/bin/filtersym ${OPT} ${FIX_FILT_OUTDIR}/${FIX_DB%.db}.${x}.${FIX_FILT_ENDING}LasFiltAln.las ${FIX_FILT_OUTDIR}/${FIX_DB%.db}.${FIX_FILT_ENDING}LasFiltAln.las.sym" 
		done > filt_09_filtsym_block_${FIX_DB%.db}.${slurmID}.plan
      	echo "DACCORD filtsym $(git --git-dir=${DACCORD_SOURCE_PATH}/.git rev-parse --short HEAD)" > filt_09_filtsym_block_${FIX_DB%.db}.${slurmID}.version                 
   	### 10_lasfilteralignmentsborderrepeats
    elif [[ ${currentStep} -eq 10 ]]
    then
        ### clean up plans 
        for x in $(ls filt_10_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        setLAfilterOptions
        
		OPT=""
        
		if [[ -z "${FIX_FILT_LASFILTERALIGNMENTSBORDERREPEATS_THREADS}" ]]
        then
        	FIX_FILT_LASFILTERALIGNMENTSBORDERREPEATS_THREADS=8
   	 	fi 
   	 	
   	 	OPT="-t${FIX_FILT_LASFILTERALIGNMENTSBORDERREPEATS_THREADS}"
   	 	
   	 	if [[ -z "${FIX_FILT_LASFILTERALIGNMENTSBORDERREPEATS_ERATE}" ]]
        then
        	FIX_FILT_LASFILTERALIGNMENTSBORDERREPEATS_ERATE=0.35
   	 	fi 
   	 	
   	 	OPT="${OPT} -e${FIX_FILT_LASFILTERALIGNMENTSBORDERREPEATS_ERATE}"
   	 	            	
    	for x in $(seq 1 ${fixblocks})
        do
    		echo "LIBMAUS2_DAZZLER_ALIGN_ALIGNMENTFILECONSTANTS_TRACE_XOVR=75 ${DACCORD_PATH}/bin/lasfilteralignmentsborderrepeats ${OPT} ${FIX_FILT_OUTDIR}/${FIX_DB%.db}.${x}.${FIX_FILT_ENDING}LasFiltBrd.las ${FIX_FILT_OUTDIR}/${FIX_DAZZ_DB%.db}.db ${FIX_FILT_OUTDIR}/${FIX_DAZZ_DB%.db}.sort.rep ${FIX_FILT_OUTDIR}/${FIX_DB%.db}.${x}.${FIX_FILT_ENDING}LasFiltAln.las" 
		done > filt_10_lasfilteralignmentsborderrepeats_block_${FIX_DB%.db}.${slurmID}.plan
      	echo "DACCORD lasfilteralignmentsborderrepeats $(git --git-dir=${DACCORD_SOURCE_PATH}/.git rev-parse --short HEAD)" > filt_10_lasfilteralignmentsborderrepeats_block_${FIX_DB%.db}.${slurmID}.version
  	### 11_mergesym2
    elif [[ ${currentStep} -eq 11 ]]
    then
        ### clean up plans 
        for x in $(ls filt_10_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        setLAfilterOptions
        OPT=""        
        echo "LIBMAUS2_DAZZLER_ALIGN_ALIGNMENTFILECONSTANTS_TRACE_XOVR=75 ${DACCORD_PATH}/bin/mergesym2 ${OPT} ${FIX_FILT_OUTDIR}/${FIX_DB%.db}.${FIX_FILT_ENDING}LasFiltBrd.las.sym ${FIX_FILT_OUTDIR}/${FIX_DAZZ_DB%.db}.db ${FIX_FILT_OUTDIR}/${FIX_DB%.db}.*.${FIX_FILT_ENDING}LasFiltBrd.las.sym" > filt_11_mergesym2_single_${FIX_DB%.db}.${slurmID}.plan
        echo "rm ${FIX_FILT_OUTDIR}/${FIX_DB%.db}.*.${FIX_FILT_ENDING}LasFiltBrd.las.sym" >> filt_11_mergesym2_single_${FIX_DB%.db}.${slurmID}.plan
        echo "DACCORD mergesym2 $(git --git-dir=${DACCORD_SOURCE_PATH}/.git rev-parse --short HEAD)" > filt_11_mergesym2_single_${FIX_DB%.db}.${slurmID}.version        
	### 12_filtersym
    elif [[ ${currentStep} -eq 12 ]]
    then
        ### clean up plans 
        for x in $(ls filt_12_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        setLAfilterOptions
        OPT=""
        
		if [[ -z "${FIX_FILT_FILTERSYM_VERBOSE}" ]]
        then
        	FIX_FILT_FILTERSYM_VERBOSE=1
   	 	fi 
   	 	
   	 	if [[ -n "${FIX_FILT_FILTERSYM_VERBOSE}" && ${FIX_FILT_FILTERSYM_VERBOSE} != 0 ]]
        then
   	 		OPT="--verbose" 
   	 	fi
   	 	
   	 	for x in $(seq 1 ${fixblocks})
        do
    		echo "LIBMAUS2_DAZZLER_ALIGN_ALIGNMENTFILECONSTANTS_TRACE_XOVR=75 ${DACCORD_PATH}/bin/filtersym ${OPT} ${FIX_FILT_OUTDIR}/${FIX_DB%.db}.${x}.${FIX_FILT_ENDING}LasFiltBrd.las ${FIX_FILT_OUTDIR}/${FIX_DB%.db}.${FIX_FILT_ENDING}LasFiltBrd.las.sym" 
		done > filt_12_filtsym_block_${FIX_DB%.db}.${slurmID}.plan
      	echo "DACCORD filtsym $(git --git-dir=${DACCORD_SOURCE_PATH}/.git rev-parse --short HEAD)" > filt_12_filtsym_block_${FIX_DB%.db}.${slurmID}.version
   	### 13_filterchainsraw
    elif [[ ${currentStep} -eq 13 ]]
    then
        ### clean up plans 
        for x in $(ls filt_13_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        setLAfilterOptions
        
        OPT=""
        
		if [[ -z "${FIX_FILT_FILTERCHAINSRAW_LEN}" ]]
        then
        	FIX_FILT_FILTERCHAINSRAW_LEN=4000
   	 	fi 
   	 	
   	 	OPT="-l${FIX_FILT_FILTERCHAINSRAW_LEN}"
        for x in $(seq 1 ${fixblocks})
        do
    		echo "LIBMAUS2_DAZZLER_ALIGN_ALIGNMENTFILECONSTANTS_TRACE_XOVR=75 ${DACCORD_PATH}/bin/filterchainsraw ${OPT} ${FIX_FILT_OUTDIR}/${FIX_DB%.db}.${x}.${FIX_FILT_ENDING}LasFiltChain.las ${FIX_FILT_OUTDIR}/${FIX_DAZZ_DB%.db}.db ${FIX_FILT_OUTDIR}/${FIX_DB%.db}.${x}.${FIX_FILT_ENDING}LasFiltBrd.las" 
		done > filt_13_filterchainsraw_block_${FIX_DB%.db}.${slurmID}.plan
        echo "DACCORD filterchainsraw $(git --git-dir=${DACCORD_SOURCE_PATH}/.git rev-parse --short HEAD)" > filt_13_filterchainsraw_block_${FIX_DB%.db}.${slurmID}.version
    ### 14_LAfilter
    elif [[ ${currentStep} -eq 14 ]]
    then
        ### clean up plans 
        for x in $(ls filt_14_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
 
        if [[ -z ${FILT_LAFILTER_OPT} ]]
        then 
            setLAfilterOptions
        fi

        if [[ -n ${FIX_FILT_LAFILTER_RMSYMROUNDS} && ${FIX_FILT_LAFILTER_RMSYMROUNDS} -gt 0 ]]
        then
            ## check what is the current round
            for rnd in $(seq ${FIX_FILT_LAFILTER_RMSYMROUNDS} -1 0)
            do
                if [[ -f filt.round${rnd}_14_LAfilter_block_${FIX_DB%.db}.${slurmID}.plan ]]
                then
                    break;
                fi
            done

            echo "stop at round $rnd of ${FIX_FILT_LAFILTER_RMSYMROUNDS}"
            ### initial filter job 
            if [[ $rnd -eq 0 ]]
            then

                ### create LAfilter commands
                for x in $(seq 1 ${fixblocks})
                do 
                    addOpt=""
                    if [[ -n ${FIX_FILT_LAFILTER_MINTIPCOV} && ${FIX_FILT_LAFILTER_MINTIPCOV} -ge 0 ]]
                    then
                        addOpt=" -a ${FIX_FILT_OUTDIR}/symDiscardOvl.round1.${x}.txt"
                    fi
                    if [[ -n ${FIX_FILT_LAFILTER_REMPERCWORSTALN} && ${FIX_FILT_LAFILTER_REMPERCWORSTALN} -gt 0 ]]
					then
    					addOpt="${addOpt} -Z ${FIX_FILT_LAFILTER_REMPERCWORSTALN}"
					fi
                    
                    echo "${MARVEL_PATH}/bin/LAfilter${FILT_LAFILTER_OPT}${addOpt} ${FIX_FILT_OUTDIR}/${FIX_DB%.db} ${FIX_FILT_OUTDIR}/${FIX_DB%.db}.${x}.${FIX_FILT_ENDING}LasFiltChain.las ${FIX_FILT_OUTDIR}/${FIX_DB%.db}.${x}.filt_R1.las"
        		done > filt.round1_14_LAfilter_block_${FIX_DB%.db}.${slurmID}.plan 
            # last filter job
            elif [[ $rnd -eq ${FIX_FILT_LAFILTER_RMSYMROUNDS} ]]
            then
                # create merged set of discarded ovls 
                cat ${FIX_FILT_OUTDIR}/symDiscardOvl.round${rnd}.*.txt | awk '{if ($1>$2) print $2" "$1; else print $1" "$2}' | sort -k1,1n -k2,2n  -u > ${FIX_FILT_OUTDIR}/symDiscardOvl.round${rnd}.txt
					
                for x in $(seq 1 ${fixblocks})
                do 
                	addOpt=""
                    if [[ -n ${FIX_FILT_LAFILTER_MINTIPCOV} && ${FIX_FILT_LAFILTER_MINTIPCOV} -ge 0 ]]
                    then
                        addOpt=" -a ${FIX_FILT_OUTDIR}/symDiscardOvl.round$((${rnd}+1)).${x}.txt -A ${FIX_FILT_OUTDIR}/symDiscardOvl.round${rnd}.txt"
                    fi
                    echo "${MARVEL_PATH}/bin/LAfilter${FILT_LAFILTER_OPT}${addOpt} ${FIX_FILT_OUTDIR}/${FIX_DB%.db} ${FIX_FILT_OUTDIR}/${FIX_DB%.db}.${x}.filt_R${rnd}.las ${FIX_FILT_OUTDIR}/${FIX_DB%.db}.${x}.filt.las"
    			done > filt_14_LAfilter_block_${FIX_DB%.db}.${slurmID}.plan 
            # intermediate filter round
            else
                # create merged set of discarded ovls 
                cat ${FIX_FILT_OUTDIR}/symDiscardOvl.round${rnd}.*.txt | awk '{if ($1>$2) print $2" "$1; else print $1" "$2}' | sort -k1,1n -k2,2n  -u > ${FIX_FILT_OUTDIR}/symDiscardOvl.round${rnd}.txt

                for x in $(seq 1 ${fixblocks})
                do 
                    addOpt=""
                    if [[ -n ${FIX_FILT_LAFILTER_MINTIPCOV} && ${FIX_FILT_LAFILTER_MINTIPCOV} -ge 0 ]]
                    then
                        addOpt=" -a ${FIX_FILT_OUTDIR}/symDiscardOvl.round$((${rnd}+1)).${x}.txt -A ${FIX_FILT_OUTDIR}/symDiscardOvl.round${rnd}.txt"
                    fi
                    echo "${MARVEL_PATH}/bin/LAfilter${FILT_LAFILTER_OPT}${addOpt} ${FIX_FILT_OUTDIR}/${FIX_DB%.db} ${FIX_FILT_OUTDIR}/${FIX_DB%.db}.${x}.filt_R${rnd}.las ${FIX_FILT_OUTDIR}/${FIX_DB%.db}.${x}.filt_R$((${rnd}+1)).las"
    		done > filt.round$((${rnd}+1))_14_LAfilter_block_${FIX_DB%.db}.${slurmID}.plan 
            fi  
        else 
            ### create LAfilter commands
            for x in $(seq 1 ${fixblocks})
            do 
                addOpt=""
                if [[ -n ${FIX_FILT_LAFILTER_MINTIPCOV} && ${FIX_FILT_LAFILTER_MINTIPCOV} -ge 0 ]]
                then
                addOpt=" -a ${FIX_FILT_OUTDIR}/discardOvlTipCov${FIX_FILT_LAFILTER_MINTIPCOV}.${x}.txt"
                fi
                echo "${MARVEL_PATH}/bin/LAfilter${FILT_LAFILTER_OPT}${addOpt} ${FIX_FILT_OUTDIR}/${FIX_DB%.db} ${FIX_FILT_OUTDIR}/${FIX_DB%.db}.${x}.${FIX_FILT_ENDING}LasFiltChain.las ${FIX_FILT_OUTDIR}/${FIX_DB%.db}.${x}.filt.las"
			done > filt_14_LAfilter_block_${FIX_DB%.db}.${slurmID}.plan 
        fi    
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > filt_14_LAfilter_block_${FIX_DB%.db}.${slurmID}.version
    #### 15_LAmerge
    elif [[ ${currentStep} -eq 15 ]]
    then
        ### clean up plans 
        for x in $(ls filt_15_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        
        if [[ -z ${FILT_LAFILTER_OPT} ]]
        then 
            setLAfilterOptions
        fi
        ### find and set LAmerge options 
        setLAmergeOptions
        
        echo "${MARVEL_PATH}/bin/LAmerge${FILT_LAMERGE_OPT} -S filt ${FIX_FILT_OUTDIR}/${FIX_DB%.db} ${FIX_FILT_OUTDIR}/${FIX_DB%.db}.filt.las" > filt_15_LAmerge_single_${FIX_DB%.db}.${slurmID}.plan
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > filt_15_LAmerge_single_${FIX_DB%.db}.${slurmID}.version
    else
        (>&2 echo "step ${currentStep} in FIX_FILT_TYPE ${FIX_FILT_TYPE} not supported")
        (>&2 echo "valid steps are: ${myTypes[${FIX_FILT_TYPE}]}")
        exit 1            
    fi
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
else
    (>&2echo "unknown RAW_PATCH_TYPE ${RAW_PATCH_TYPE}")    
    (>&2 echo "supported types")
    x=0; while [ $x -lt ${#myTypes[*]} ]; do (>&2 echo "${myTypes[${x}]}"); done
    exit 1
fi

exit 0
