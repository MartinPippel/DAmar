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

function setLAfilterOptions()
{
	FIX_LAFILTER_OPT=""
	
	if [[ -z "${RAW_FIX_LAFILTER_PURGE}" ]]
	then
		RAW_FIX_LAFILTER_PURGE=1	
	fi 
	
	if [[ ${RAW_FIX_LAFILTERCHAINS_ANCHOR} -gt 0 ]]
	then 
		FIX_LAFILTER_OPT="${FIX_LAFILTER_OPT} -p"
	fi
	
	if [[ -z "${RAW_FIX_LAFILTER_MAXSEGERR}" ]]
	then
		RAW_FIX_LAFILTER_MAXSEGERR=65	
	fi 
	
	if [[ ${RAW_FIX_LAFILTER_MAXSEGERR} -gt 0 ]]
	then 
		FIX_LAFILTER_OPT="${FIX_LAFILTER_OPT} -b ${RAW_FIX_LAFILTER_MAXSEGERR}"
	fi
}

function setLAfilterChainsOptions()
{
	FIX_LAFILTERCHAINS_OPT=""
	
	if [[ -n ${RAW_FIX_LAFILTERCHAINS_ANCHOR} && ${RAW_FIX_LAFILTERCHAINS_ANCHOR} -gt 0 ]]
	then 
		FIX_LAFILTERCHAINS_OPT="${FIX_LAFILTERCHAINS_OPT} -n ${RAW_FIX_LAFILTERCHAINS_ANCHOR}"
	fi
	
	if [[ -n ${RAW_FIX_LAFILTERCHAINS_PURGE} && ${RAW_FIX_LAFILTERCHAINS_PURGE} -gt 0 ]]
	then 
		FIX_LAFILTERCHAINS_OPT="${FIX_LAFILTERCHAINS_OPT} -p"
	fi
	
	if [[ -n ${RAW_FIX_LAFILTERCHAINS_NKEEPCHAINS} ]]
	then 
		FIX_LAFILTERCHAINS_OPT="${FIX_LAFILTERCHAINS_OPT} -k ${RAW_FIX_LAFILTERCHAINS_NKEEPCHAINS}"
	fi
	
	if [[ -n ${RAW_FIX_LAFILTERCHAINS_LOWCOMP} ]]
	then 
		FIX_LAFILTERCHAINS_OPT="${FIX_LAFILTERCHAINS_OPT} -l ${RAW_FIX_LAFILTERCHAINS_LOWCOMP}"
	fi
	
	# add default trim and q tracks
	
	if [[ "${RAW_DACCORD_INDIR}" == "${RAW_REPCOMP_OUTDIR}" ]]
	then
		FIX_LAFILTERCHAINS_OPT="${FIX_LAFILTERCHAINS_OPT} -t trim0_d${RAW_FIX_LAQ_QTRIMCUTOFF}_s${RAW_FIX_LAQ_MINSEG}_repcomp -q q0_d${RAW_FIX_LAQ_QTRIMCUTOFF}_s${RAW_FIX_LAQ_MINSEG}_repcomp"	
	else
		FIX_LAFILTERCHAINS_OPT="${FIX_LAFILTERCHAINS_OPT} -t trim0_d${RAW_FIX_LAQ_QTRIMCUTOFF}_s${RAW_FIX_LAQ_MINSEG}_dalign -q q0_d${RAW_FIX_LAQ_QTRIMCUTOFF}_s${RAW_FIX_LAQ_MINSEG}_dalign"
	fi
	
	if [[ -n ${RAW_FIX_LAFILTERCHAINS_UNALIGNBASES} && ${RAW_FIX_LAFILTERCHAINS_UNALIGNBASES} -gt 0 ]]
	then 
		FIX_LAFILTERCHAINS_OPT="${FIX_LAFILTERCHAINS_OPT} -u ${RAW_FIX_LAFILTERCHAINS_UNALIGNBASES}"
	fi
	
	if [[ -n ${RAW_FIX_LAFILTERCHAINS_DIFF} && ${RAW_FIX_LAFILTERCHAINS_DIFF} -gt 0 ]]
	then 
		FIX_LAFILTERCHAINS_OPT="${FIX_LAFILTERCHAINS_OPT} -d ${RAW_FIX_LAFILTERCHAINS_DIFF}"
	fi
	
	if [[ -n ${RAW_FIX_LAFILTERCHAINS_CHAINLEN} && ${RAW_FIX_LAFILTERCHAINS_CHAINLEN} -gt 0 ]]
	then 
		FIX_LAFILTERCHAINS_OPT="${FIX_LAFILTERCHAINS_OPT} -o ${RAW_FIX_LAFILTERCHAINS_CHAINLEN}"
	fi
	
	if [[ -n ${RAW_FIX_LAFILTERCHAINS_FULLDISCAREADLEN} && ${RAW_FIX_LAFILTERCHAINS_FULLDISCAREADLEN} -gt 0 ]]
	then 
		FIX_LAFILTERCHAINS_OPT="${FIX_LAFILTERCHAINS_OPT} -Z ${RAW_FIX_LAFILTERCHAINS_FULLDISCAREADLEN}"
	fi
	
}

function setDaccordOptions()
{
	FIX_DACCORD_OPT=""
	
	if [[ -z ${RAW_FIX_DACCORD_THREADS} ]]
	then 
		RAW_FIX_DACCORD_THREADS=8
	fi
	FIX_DACCORD_OPT="${FIX_DACCORD_OPT} -t${RAW_FIX_DACCORD_THREADS}"
	
	if [[ -n ${RAW_FIX_DACCORD_WINDOW} && ${RAW_FIX_DACCORD_WINDOW} -gt 0 ]]
	then 
		FIX_DACCORD_OPT="${FIX_DACCORD_OPT} -w${RAW_FIX_DACCORD_WINDOW}"
	fi

	if [[ -n ${RAW_FIX_DACCORD_ADVANCESIZE} && ${RAW_FIX_DACCORD_ADVANCESIZE} -gt 0 ]]
	then 
		FIX_DACCORD_OPT="${FIX_DACCORD_OPT} -a${RAW_FIX_DACCORD_ADVANCESIZE}"
	fi
	
	if [[ -n ${RAW_FIX_DACCORD_MAXDEPTH} && ${RAW_FIX_DACCORD_MAXDEPTH} -gt 0 ]]
	then 
		FIX_DACCORD_OPT="${FIX_DACCORD_OPT} -d${RAW_FIX_DACCORD_MAXDEPTH}"
	fi
	
	if [[ -n ${RAW_FIX_DACCORD_FULLSEQ} && ${RAW_FIX_DACCORD_FULLSEQ} -gt 0 ]]
	then 
		FIX_DACCORD_OPT="${FIX_DACCORD_OPT} -f1"
	fi
	
	if [[ -n ${RAW_FIX_DACCORD_VEBOSE} && ${RAW_FIX_DACCORD_VEBOSE} -gt 0 ]]
	then 
		FIX_DACCORD_OPT="${FIX_DACCORD_OPT} -V${RAW_FIX_DACCORD_VEBOSE}"
	fi
		
	if [[ -n ${RAW_FIX_DACCORD_MINWINDOWCOV} && ${RAW_FIX_DACCORD_MINWINDOWCOV} -gt 0 ]]
	then 
		FIX_DACCORD_OPT="${FIX_DACCORD_OPT} -m${RAW_FIX_DACCORD_MINWINDOWCOV}"
	fi
	
	if [[ -n ${RAW_FIX_DACCORD_MINWINDOWERR} && ${RAW_FIX_DACCORD_MINWINDOWERR} -gt 0 ]]
	then 
		FIX_DACCORD_OPT="${FIX_DACCORD_OPT} -e${RAW_FIX_DACCORD_MINWINDOWERR}"
	fi
	
	if [[ -n ${RAW_FIX_DACCORD_MINOUTLEN} && ${RAW_FIX_DACCORD_MINOUTLEN} -gt 0 ]]
	then 
		FIX_DACCORD_OPT="${FIX_DACCORD_OPT} -l${RAW_FIX_DACCORD_MINOUTLEN}"
	fi
	
	if [[ -n ${RAW_FIX_DACCORD_MINKFREQ} && ${RAW_FIX_DACCORD_MINKFREQ} -gt 0 ]]
	then 
		FIX_DACCORD_OPT="${FIX_DACCORD_OPT} --minfilterfreq${RAW_FIX_DACCORD_MINKFREQ}"
	fi
	
	if [[ -n ${RAW_FIX_DACCORD_MAXKFREQ} && ${RAW_FIX_DACCORD_MAXKFREQ} -gt 0 ]]
	then 
		FIX_DACCORD_OPT="${FIX_DACCORD_OPT} --maxfilterfreq${RAW_FIX_DACCORD_MAXKFREQ}"
	fi
	
	if [[ -n ${RAW_FIX_DACCORD_MAXOVLS} && ${RAW_FIX_DACCORD_MAXOVLS} -gt 0 ]]
	then 
		FIX_DACCORD_OPT="${FIX_DACCORD_OPT} -D${RAW_FIX_DACCORD_MAXOVLS}"
	fi
	
	if [[ -n ${RAW_FIX_DACCORD_VARD} && ${RAW_FIX_DACCORD_VARD} -gt 0 ]]
	then 
		FIX_DACCORD_OPT="${FIX_DACCORD_OPT} --vard${RAW_FIX_DACCORD_VARD}"
	fi
	
	if [[ -n ${RAW_FIX_DACCORD_KMER} && ${RAW_FIX_DACCORD_KMER} -gt 0 ]]
	then 
		FIX_DACCORD_OPT="${FIX_DACCORD_OPT} -k${RAW_FIX_DACCORD_KMER}"
	fi
}

function setHaploSplitOptions()
{
	FIX_SPLIT_OPT=""
	
	if [[ -z "${RAW_FIX_SPLIT_TYPE}" ]]
	then 
		RAW_FIX_SPLIT_TYPE=split_agr		
	fi
	
	if [[ "${RAW_FIX_SPLIT_TYPE}" != "split_agr" && "${RAW_FIX_SPLIT_TYPE}" != "split_dis" ]]
	then
		(>&2 echo "ERROR - Split type not supported. (split_agr or split_dis)")
    	exit 1
	fi
	
	FIX_SPLIT_OPT="${RAW_FIX_SPLIT_TYPE}"
	
	if [[ -z "${RAW_FIX_SPLIT_THREADS}" ]]
	then
		RAW_FIX_SPLIT_THREADS=8
	fi
	
	FIX_SPLIT_OPT="${FIX_SPLIT_OPT} -t${RAW_FIX_SPLIT_THREADS}"	
	
	if [[ -n "${RAW_FIX_SPLIT_SEQDEPTH}" && ${RAW_FIX_SPLIT_SEQDEPTH} -gt 0 ]]
	then
		FIX_SPLIT_OPT="${FIX_SPLIT_OPT} -d${RAW_FIX_SPLIT_SEQDEPTH}"
	fi
	
	if [[ -n "${RAW_FIX_SPLIT_PHASETHRESHOLD}" ]]
	then
		FIX_SPLIT_OPT="${FIX_SPLIT_OPT} -p${RAW_FIX_SPLIT_PHASETHRESHOLD}"
	fi
	
	if [[ -n "${RAW_FIX_SPLIT_MAXALNS}" && ${RAW_FIX_SPLIT_MAXALNS} -gt 0 ]]
	then
		FIX_SPLIT_OPT="${FIX_SPLIT_OPT} -D${RAW_FIX_SPLIT_MAXALNS}"
	fi
		
	### for split split_dis there are further options available
	if [[ -n ${RAW_FIX_SPLIT_DIFFRATE} ]]
	then 
		FIX_SPLIT_OPT="${FIX_SPLIT_OPT} --drate${RAW_FIX_SPLIT_DIFFRATE}"
	fi
	
	if [[ -n ${RAW_FIX_SPLIT_NUMVARS} ]]
	then 
		FIX_SPLIT_OPT="${FIX_SPLIT_OPT} --kv${RAW_FIX_SPLIT_NUMVARS}"
	fi

	if [[ -n ${RAW_FIX_SPLIT_PHASETYPE} ]]
	then 
		FIX_SPLIT_OPT="${FIX_SPLIT_OPT} --phasetype${RAW_FIX_SPLIT_PHASETYPE}"
	fi
	
	
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

if [[ ${currentStep} -lt 10 ]]
then 
	sID=0${currentStep}
else
	sID=${currentStep}
fi
myCWD=$(pwd)

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
        for x in $(ls fix_${sID}_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        
        echo "if [[ -d ${RAW_DALIGN_OUTDIR} ]]; then mv ${RAW_DALIGN_OUTDIR} ${RAW_DALIGN_OUTDIR}_$(date '+%Y-%m-%d_%H-%M-%S'); fi && mkdir ${RAW_DALIGN_OUTDIR} && ln -s -r .${RAW_DB%.db}.* ${RAW_DB%.db}.db .${RAW_DAZZ_DB%.db}.* ${RAW_DAZZ_DB%.db}.db ${RAW_DALIGN_OUTDIR}" > fix_${sID}_createSubdir_single_${RAW_DB%.db}.${slurmID}.plan
		for x in $(seq 1 ${nblocks})
	        do
			echo "mkdir -p ${RAW_DALIGN_OUTDIR}/d${x}"
		done >> fix_${sID}_createSubdir_single_${RAW_DB%.db}.${slurmID}.plan
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_${sID}_createSubdir_single_${RAW_DB%.db}.${slurmID}.version
    elif [[ ${currentStep} -eq 2 ]]
    then
        ### clean up plans 
        for x in $(ls fix_${sID}_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        
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
            	if [[ "x${DALIGNER_VERSION}" == "x2" ]]
            	then
            		echo -n "cd ${RAW_DALIGN_OUTDIR} && PATH=${DAZZLER_PATH}/bin:\${PATH} ${NUMACTL}${DAZZLER_PATH}/bin/daligner${FIX_DALIGNER_OPT} ${RAW_DAZZ_DB%.db}.${x} ${RAW_DAZZ_DB%.db}.@${x}"
		else
        		echo -n "cd ${RAW_DALIGN_OUTDIR} && PATH=${DAZZLER_PATH}/bin:\${PATH} ${NUMACTL}${DAZZLER_PATH}/bin/daligner${FIX_DALIGNER_OPT} ${RAW_DAZZ_DB%.db}.${x}"
		fi
            cmdLine=$((${cmdLine}+1))
            count=0

            for y in $(seq ${x} ${nblocks})
            do  
                if [[ $count -lt ${RAW_FIX_DALIGNER_DAL} ]]
                then
                    count=$((${count}+1))
                    echo -n " ${RAW_DAZZ_DB%.db}.${y}"
                else
                	if [[ "x${DALIGNER_VERSION}" == "x2" ]]
            		then    
                    	        echo -n "-$((y-1)) && mv"
                	else
                		echo -n " && mv"
                	fi
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
                    if [[ "x${DALIGNER_VERSION}" == "x2" ]]
            		then
                    		echo -n "cd ${RAW_DALIGN_OUTDIR} && PATH=${DAZZLER_PATH}/bin:\${PATH} ${NUMACTL}${DAZZLER_PATH}/bin/daligner${FIX_DALIGNER_OPT} ${RAW_DAZZ_DB%.db}.${x} ${RAW_DAZZ_DB%.db}.@${y}"
                	else
                		echo -n "cd ${RAW_DALIGN_OUTDIR} && PATH=${DAZZLER_PATH}/bin:\${PATH} ${NUMACTL}${DAZZLER_PATH}/bin/daligner${FIX_DALIGNER_OPT} ${RAW_DAZZ_DB%.db}.${x} ${RAW_DAZZ_DB%.db}.${y}"
                	fi
                    cmdLine=$((${cmdLine}+1))
                    count=1
                fi
            done
	    if [[ "x${DALIGNER_VERSION}" == "x2" ]]	
	    then
            	echo -n "-${y} && mv"
	    else
		echo -n " && mv"
	    fi
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
    	done > fix_${sID}_daligner_block_${RAW_DB%.db}.${slurmID}.plan
        echo "DAZZLER daligner $(git --git-dir=${DAZZLER_SOURCE_PATH}/DALIGNER/.git rev-parse --short HEAD)" > fix_${sID}_daligner_block_${RAW_DB%.db}.${slurmID}.version
    elif [[ ${currentStep} -eq 3 ]]
    then
        ### clean up plans 
        for x in $(ls fix_${sID}_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        ### find and set LAmerge options 
        setLAmergeOptions
        ### create LAmerge commands
        for x in $(seq 1 ${nblocks})
        do 
            echo "cd ${RAW_DALIGN_OUTDIR} && ${MARVEL_PATH}/bin/LAmerge${FIX_LAMERGE_OPT} ${RAW_DB%.db} ${RAW_DAZZ_DB%.db}.dalign.${x}.las d${x} && ${MARVEL_PATH}/bin/LAfilter -p -R6 ${RAW_DB%.db} ${RAW_DAZZ_DB%.db}.dalign.${x}.las ${RAW_DAZZ_DB%.db}.dalignFilt.${x}.las && cd ${myCWD}"
    	done > fix_${sID}_LAmerge_block_${RAW_DB%.db}.${slurmID}.plan  
        echo "MARVEL LAmerge $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_${sID}_LAmerge_block_${RAW_DB%.db}.${slurmID}.version       
    elif [[ ${currentStep} -eq 4 ]]
    then
        ### clean up plans 
        for x in $(ls fix_${sID}_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        setLArepeatOptions 1
        ### create LArepeat commands
        for x in $(seq 1 ${nblocks})
        do 
            echo "cd ${RAW_DALIGN_OUTDIR} && ${MARVEL_PATH}/bin/LArepeat${FIX_LAREPEAT_OPT} -b ${x} ${RAW_DB%.db} ${RAW_DAZZ_DB%.db}.dalignFilt.${x}.las && cd ${myCWD}/"
            echo "cd ${RAW_DALIGN_OUTDIR} && ${DAZZLER_PATH}/bin/REPmask -v -c${RAW_DAZZ_FIX_LAREPEAT_THRESHOLD} -n${RAW_DAZZ_FIX_LAREPEAT_REPEATTRACK} ${RAW_DAZZ_DB%.db} ${RAW_DAZZ_DB%.db}.dalignFilt.${x}.las && cd ${myCWD}/"
    	done > fix_${sID}_LArepeat_block_${RAW_DB%.db}.${slurmID}.plan 
        echo "MARVEL LArepeat $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_${sID}_LArepeat_block_${RAW_DB%.db}.${slurmID}.version
        echo "DAZZLER REPmask $(git --git-dir=${DAZZLER_SOURCE_PATH}/DAMASKER/.git rev-parse --short HEAD)" >> fix_${sID}_LArepeat_block_${RAW_DB%.db}.${slurmID}.version        
    elif [[ ${currentStep} -eq 5 ]]
    then
        ### clean up plans 
        for x in $(ls fix_${sID}_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        
		# we need the name of the repeat track, especially if the plan starts with step4
        setLArepeatOptions 1
        ### find and set TKmerge options 
        if [[ -z "${FIX_TKMERGE_OPT}" ]]
        then 
            setTKmergeOptions
        fi
               
        ### create TKmerge command
        echo "cd ${RAW_DALIGN_OUTDIR} && ${MARVEL_PATH}/bin/TKmerge${FIX_TKMERGE_OPT} ${RAW_DB%.db} ${RAW_FIX_LAREPEAT_REPEATTRACK} && cp .${RAW_DB%.db}.${RAW_FIX_LAREPEAT_REPEATTRACK}.a2 .${RAW_DB%.db}.${RAW_FIX_LAREPEAT_REPEATTRACK}.d2 ${myCWD} && cd ${myCWD}" > fix_${sID}_TKmerge_single_${RAW_DB%.db}.${slurmID}.plan      
        echo "cd ${RAW_DALIGN_OUTDIR} && ${DAZZLER_PATH}/bin/Catrack${FIX_TKMERGE_OPT} -f -v ${RAW_DAZZ_DB%.db} ${RAW_DAZZ_FIX_LAREPEAT_REPEATTRACK} && cp .${RAW_DAZZ_DB%.db}.${RAW_DAZZ_FIX_LAREPEAT_REPEATTRACK}.anno .${RAW_DAZZ_DB%.db}.${RAW_DAZZ_FIX_LAREPEAT_REPEATTRACK}.data ${myCWD}/ && cd ${myCWD}/" >> fix_${sID}_TKmerge_single_${RAW_DB%.db}.${slurmID}.plan
        echo "MARVEL TKmerge $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_${sID}_TKmerge_single_${RAW_DB%.db}.${slurmID}.version
        echo "DAZZLER Catrack $(git --git-dir=${DAZZLER_SOURCE_PATH}/DAZZ_DB/.git rev-parse --short HEAD)" >> fix_${sID}_TKmerge_single_${RAW_DB%.db}.${slurmID}.version   
    elif [[ ${currentStep} -eq 6 ]]
    then
        ### clean up plans 
        for x in $(ls fix_${sID}_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
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
            echo "${MARVEL_PATH}/bin/TKcombine${FIX_TKCOMBINE_OPT} ${RAW_DB%.db} ${RAW_FIX_LAREPEAT_REPEATTRACK}_${RAW_REPMASK_LAREPEAT_REPEATTRACK} ${RAW_FIX_LAREPEAT_REPEATTRACK} ${RAW_REPMASK_REPEATTRACK}" > fix_${sID}_TKcombine_single_${RAW_DB%.db}.${slurmID}.plan
            echo "${MARVEL_PATH}/bin/TKcombine${FIX_TKCOMBINE_OPT} ${RAW_DB%.db} ${RAW_FIX_LAREPEAT_REPEATTRACK}_${RAW_REPMASK_LAREPEAT_REPEATTRACK}_${RAW_REPMASK_TANMASK_TRACK}_dust ${RAW_FIX_LAREPEAT_REPEATTRACK}_${RAW_REPMASK_LAREPEAT_REPEATTRACK} ${RAW_REPMASK_TANMASK_TRACK} dust" >> fix_${sID}_TKcombine_single_${RAW_DB%.db}.${slurmID}.plan         
        else
            echo "ln -s .${RAW_DB%.db}.${RAW_FIX_LAREPEAT_REPEATTRACK}.d2 .${RAW_DB%.db}.${RAW_FIX_LAREPEAT_REPEATTRACK}_${RAW_REPMASK_LAREPEAT_REPEATTRACK}.d2"  > fix_${sID}_TKcombine_single_${RAW_DB%.db}.${slurmID}.plan         
            echo "ln -s .${RAW_DB%.db}.${RAW_FIX_LAREPEAT_REPEATTRACK}.a2 .${RAW_DB%.db}.${RAW_FIX_LAREPEAT_REPEATTRACK}_${RAW_REPMASK_LAREPEAT_REPEATTRACK}.a2"  >> fix_${sID}_TKcombine_single_${RAW_DB%.db}.${slurmID}.plan         
        fi 
        echo "MARVEL TKcombine $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_${sID}_TKcombine_single_${RAW_DB%.db}.${slurmID}.version         
    elif [[ ${currentStep} -eq 7 ]]
    then
        ### clean up plans 
        for x in $(ls fix_${sID}_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 

        mkdir -p identity
        ### create LAfilter commands - filter out identity overlaps - has to be done because revcomp and forcealign will loose those 
        for x in $(seq 1 ${nblocks})
        do  
            echo "${MARVEL_PATH}/bin/LAfilter -p -R 3 -R 6 ${RAW_DB%.db} ${RAW_DALIGN_OUTDIR}/d${x}/${RAW_DAZZ_DB%.db}.${x}.${RAW_DAZZ_DB%.db}.${x}.las identity/${RAW_DAZZ_DB%.db}.identity.${x}.las"
		done > fix_${sID}_LAfilter_block_${RAW_DB%.db}.${slurmID}.plan   
    echo "MARVEL LAfilter $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_${sID}_LAfilter_block_${RAW_DB%.db}.${slurmID}.version    
    elif [[ ${currentStep} -eq 8 ]]
    then
        ### clean up plans 
        for x in $(ls fix_${sID}_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        
		### find and set LAq options 
        setLAqOptions
        ### create LAq commands
        for x in $(seq 1 ${nblocks})
        do 
            echo "cd ${RAW_DALIGN_OUTDIR} && ${MARVEL_PATH}/bin/LAq${FIX_LAQ_OPT} -T trim0_d${RAW_FIX_LAQ_QTRIMCUTOFF}_s${RAW_FIX_LAQ_MINSEG}_dalign -Q q0_d${RAW_FIX_LAQ_QTRIMCUTOFF}_s${RAW_FIX_LAQ_MINSEG}_dalign ${RAW_DB%.db} -b ${x} ${RAW_DAZZ_DB%.db}.dalignFilt.${x}.las && cd ${myCWD}"
    	done > fix_${sID}_LAq_block_${RAW_DB%.db}.${slurmID}.plan 
        echo "MARVEL LAq $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_${sID}_LAq_block_${RAW_DB%.db}.${slurmID}.version                
    elif [[ ${currentStep} -eq 9 ]]
    then
        ### clean up plans 
        for x in $(ls fix_${sID}_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
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
        
		### create TKmerge command
        echo "cd ${RAW_DALIGN_OUTDIR} && ${MARVEL_PATH}/bin/TKmerge${FIX_TKMERGE_OPT} ${RAW_DB%.db} trim0_d${RAW_FIX_LAQ_QTRIMCUTOFF}_s${RAW_FIX_LAQ_MINSEG}_dalign && cp .${RAW_DB%.db}.trim0_d${RAW_FIX_LAQ_QTRIMCUTOFF}_s${RAW_FIX_LAQ_MINSEG}_dalign.a2 .${RAW_DB%.db}.trim0_d${RAW_FIX_LAQ_QTRIMCUTOFF}_s${RAW_FIX_LAQ_MINSEG}_dalign.d2 ${myCWD}/ && cd ${myCWD}" > fix_${sID}_TKmerge_block_${RAW_DB%.db}.${slurmID}.plan
        echo "cd ${RAW_DALIGN_OUTDIR} && ${MARVEL_PATH}/bin/TKmerge${FIX_TKMERGE_OPT} ${RAW_DB%.db} q0_d${RAW_FIX_LAQ_QTRIMCUTOFF}_s${RAW_FIX_LAQ_MINSEG}_dalign && cp .${RAW_DB%.db}.q0_d${RAW_FIX_LAQ_QTRIMCUTOFF}_s${RAW_FIX_LAQ_MINSEG}_dalign.a2 .${RAW_DB%.db}.q0_d${RAW_FIX_LAQ_QTRIMCUTOFF}_s${RAW_FIX_LAQ_MINSEG}_dalign.d2 ${myCWD}/ && cd ${myCWD}" >> fix_${sID}_TKmerge_block_${RAW_DB%.db}.${slurmID}.plan       
        echo "MARVEL TKmerge $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > ${sID}_TKmerge_block_${RAW_DB%.db}.${slurmID}.version               
	### 10_LAfix    
    elif [[ ${currentStep} -eq 10 ]]
    then
        ### clean up plans 
        for x in $(ls fix_${sID}_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
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
            echo "${MARVEL_PATH}/bin/LAfix${FIX_LAFIX_OPT} ${addopt}${RAW_DB%.db} ${RAW_DALIGN_OUTDIR}/${RAW_DAZZ_DB%.db}.dalignFilt.${x}.las ${RAW_FIX_LAFIX_PATH}/${RAW_DB%.db}.${x}${RAW_FIX_LAFIX_FILESUFFIX}.fasta"
    	done > fix_${sID}_LAfix_block_${RAW_DB%.db}.${slurmID}.plan
    echo "MARVEL LAfix $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_${sID}_LAfix_block_${RAW_DB%.db}.${slurmID}.version                 
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
	    for x in $(ls fix_${sID}_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
	    do            
	        rm $x
	    done 
	    
	    echo "if [[ -d ${RAW_REPCOMP_OUTDIR} ]]; then mv ${RAW_REPCOMP_OUTDIR} ${RAW_REPCOMP_OUTDIR}_$(date '+%Y-%m-%d_%H-%M-%S'); fi && mkdir ${RAW_REPCOMP_OUTDIR} && ln -s -r .${RAW_DB%.db}.* ${RAW_DB%.db}.db .${RAW_DAZZ_DB%.db}.* ${RAW_DAZZ_DB%.db}.db ${RAW_REPCOMP_OUTDIR}" > fix_${sID}_createSubdir_single_${RAW_DB%.db}.${slurmID}.plan
		for x in $(seq 1 ${nblocks})
	    do
			echo "mkdir -p ${RAW_REPCOMP_OUTDIR}/r${x} ${RAW_REPCOMP_OUTDIR}/d${x}_ForRepComp ${RAW_REPCOMP_OUTDIR}/d${x}_NoRepComp"
		done >> fix_${sID}_createSubdir_single_${RAW_DB%.db}.${slurmID}.plan
	    echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_${sID}_createSubdir_single_${RAW_DB%.db}.${slurmID}.version
  	#### LAseparate
  	elif [[ ${currentStep} -eq 2 ]]
    then
        ### clean up plans 
        for x in $(ls fix_${sID}_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
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
                echo "${MARVEL_PATH}/bin/LAseparate${FIX_LASEPARATE_OPT} ${RAW_DB%.db} ${RAW_DALIGN_OUTDIR}/d${x}/${RAW_DAZZ_DB%.db}.${x}.${RAW_DAZZ_DB%.db}.${y}.las ${RAW_REPCOMP_OUTDIR}/d${x}_ForRepComp/${RAW_DAZZ_DB%.db}.${x}.${RAW_DAZZ_DB%.db}.${y}.las ${RAW_REPCOMP_OUTDIR}/d${x}_NoRepComp/${RAW_DAZZ_DB%.db}.${x}.${RAW_DAZZ_DB%.db}.${y}.las"                
            done 
    	done > fix_${sID}_LAseparate_block_${RAW_DB%.db}.${slurmID}.plan
    	echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_${sID}_LAseparate_block_${RAW_DB%.db}.${slurmID}.version
    #### repcomp 
    elif [[ ${currentStep} -eq 3 ]]
    then
        ### clean up plans 
        for x in $(ls fix_${sID}_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
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
		done > fix_${sID}_repcomp_block_${RAW_DB%.db}.${slurmID}.plan
    	echo "repcomp $(git --git-dir=${REPCOMP_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_${sID}_repcomp_block_${RAW_DB%.db}.${slurmID}.version
	### 04_LAmergeLAfilter
    elif [[ ${currentStep} -eq 4 ]]
    then
        ### clean up plans 
        for x in $(ls fix_${sID}_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        
        ### find and set LAmerge options 
        setLAmergeOptions
        setRepcompOptions
        ### create LAmerge commands
        for x in $(seq 1 ${nblocks})
        do 
            echo "cd ${RAW_REPCOMP_OUTDIR} && ${MARVEL_PATH}/bin/LAmerge${FIX_LAMERGE_OPT} ${RAW_DB%.db} ${RAW_DAZZ_DB%.db}.repcomp.${x}.las r${x} d${x}_ForRepComp d${x}_NoRepComp ${myCWD}/identity/${RAW_DAZZ_DB%.db}.identity.${x}.las && ${MARVEL_PATH}/bin/LAfilter -p -R6 ${RAW_DB%.db} ${RAW_DAZZ_DB%.db}.repcomp.${x}.las ${RAW_DAZZ_DB%.db}.repcompFilt.${x}.las && cd ${myCWD}"                                                                                                                     
    	done > fix_${sID}_LAmerge_block_${RAW_DB%.db}.${slurmID}.plan
    	echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_${sID}_LAmerge_block_${RAW_DB%.db}.${slurmID}.version      
    ### 05_LArepeat
    elif [[ ${currentStep} -eq 5 ]]
    then
        ### clean up plans 
        for x in $(ls fix_${sID}_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        
        setLArepeatOptions 2
        ### create LArepeat commands
        for x in $(seq 1 ${nblocks})
        do 
            echo "cd ${RAW_REPCOMP_OUTDIR} && ${MARVEL_PATH}/bin/LArepeat${FIX_LAREPEAT_OPT} -b ${x} ${RAW_DB%.db} ${RAW_DAZZ_DB%.db}.repcompFilt.${x}.las && cd ${myCWD}/"
            echo "cd ${RAW_REPCOMP_OUTDIR} && ${DAZZLER_PATH}/bin/REPmask -v -c${RAW_DAZZ_FIX_LAREPEAT_THRESHOLD} -n${RAW_DAZZ_FIX_LAREPEAT_REPEATTRACK} ${RAW_DAZZ_DB%.db} ${RAW_DAZZ_DB%.db}.repcompFilt.${x}.las && cd ${myCWD}/"
            
    	done > fix_${sID}_LArepeat_block_${RAW_DB%.db}.${slurmID}.plan
    	echo "MARVEL LArepeat $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_${sID}_LArepeat_block_${RAW_DB%.db}.${slurmID}.version
    	echo "DAZZLER REPmask $(git --git-dir=${DAZZLER_SOURCE_PATH}/DAMASKER/.git rev-parse --short HEAD)" >> fix_${sID}_LArepeat_block_${RAW_DB%.db}.${slurmID}.version
    ### 06_TKmerge         
    elif [[ ${currentStep} -eq 6 ]]
    then
        ### clean up plans 
        for x in $(ls fix_${sID}_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        
        # we need the name of the repeat track, especially if the plan starts with step4
        setLArepeatOptions 2
        ### find and set TKmerge options 
        if [[ -z ${FIX_TKMERGE_OPT} ]]
        then 
            setTKmergeOptions
        fi
        ### create TKmerge command
        echo "cd ${RAW_REPCOMP_OUTDIR} && ${MARVEL_PATH}/bin/TKmerge${FIX_TKMERGE_OPT} ${RAW_DB%.db} ${RAW_FIX_LAREPEAT_REPEATTRACK} && cp .${RAW_DB%.db}.${RAW_FIX_LAREPEAT_REPEATTRACK}.a2 .${RAW_DB%.db}.${RAW_FIX_LAREPEAT_REPEATTRACK}.d2 ${myCWD} && cd ${myCWD}" > fix_${sID}_TKmerge_single_${RAW_DB%.db}.${slurmID}.plan      
        echo "cd ${RAW_REPCOMP_OUTDIR} && ${DAZZLER_PATH}/bin/Catrack${FIX_TKMERGE_OPT} -f -v ${RAW_DAZZ_DB%.db} ${RAW_DAZZ_FIX_LAREPEAT_REPEATTRACK} && cp .${RAW_DAZZ_DB%.db}.${RAW_DAZZ_FIX_LAREPEAT_REPEATTRACK}.anno .${RAW_DAZZ_DB%.db}.${RAW_DAZZ_FIX_LAREPEAT_REPEATTRACK}.data ${myCWD}/ && cd ${myCWD}/" >> fix_${sID}_TKmerge_single_${RAW_DB%.db}.${slurmID}.plan
        echo "MARVEL TKmerge $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_${sID}_TKmerge_single_${RAW_DB%.db}.${slurmID}.version
        echo "DAZZLER Catrack $(git --git-dir=${DAZZLER_SOURCE_PATH}/DAZZ_DB/.git rev-parse --short HEAD)" >> fix_${sID}_TKmerge_single_${RAW_DB%.db}.${slurmID}.version
    ### 07_TKcombine   
    elif [[ ${currentStep} -eq 7 ]]
    then
        ### clean up plans 
        for x in $(ls fix_${sID}_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
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
            echo "${MARVEL_PATH}/bin/TKcombine${FIX_TKCOMBINE_OPT} ${RAW_DB%.db} ${RAW_FIX_LAREPEAT_REPEATTRACK}_${RAW_REPMASK_LAREPEAT_REPEATTRACK} ${RAW_FIX_LAREPEAT_REPEATTRACK} ${RAW_REPMASK_REPEATTRACK}" > fix_${sID}_TKcombine_single_${RAW_DB%.db}.${slurmID}.plan
            echo "${MARVEL_PATH}/bin/TKcombine${FIX_TKCOMBINE_OPT} ${RAW_DB%.db} ${RAW_FIX_LAREPEAT_REPEATTRACK}_${RAW_REPMASK_LAREPEAT_REPEATTRACK}_${RAW_REPMASK_TANMASK_TRACK}_dust ${RAW_FIX_LAREPEAT_REPEATTRACK}_${RAW_REPMASK_LAREPEAT_REPEATTRACK} ${RAW_REPMASK_TANMASK_TRACK} dust" >> fix_${sID}_TKcombine_single_${RAW_DB%.db}.${slurmID}.plan         
        else
            echo "ln -s .${RAW_DB%.db}.${RAW_FIX_LAREPEAT_REPEATTRACK}.d2 .${RAW_DB%.db}.${RAW_FIX_LAREPEAT_REPEATTRACK}_${RAW_REPMASK_LAREPEAT_REPEATTRACK}.d2"  > fix_${sID}_TKcombine_single_${RAW_DB%.db}.${slurmID}.plan         
            echo "ln -s .${RAW_DB%.db}.${RAW_FIX_LAREPEAT_REPEATTRACK}.a2 .${RAW_DB%.db}.${RAW_FIX_LAREPEAT_REPEATTRACK}_${RAW_REPMASK_LAREPEAT_REPEATTRACK}.a2"  >> fix_${sID}_TKcombine_single_${RAW_DB%.db}.${slurmID}.plan         
        fi 
        echo "MARVEL TKcombine $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_${sID}_TKcombine_single_${RAW_DB%.db}.${slurmID}.version
    ### 08_LAq  
    elif [[ ${currentStep} -eq 8 ]]
    then
        ### clean up plans 
        for x in $(ls fix_${sID}_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        
        ### find and set LAq options 
        setLAqOptions
        ### create LAq commands
        for x in $(seq 1 ${nblocks})
        do 
            echo "cd ${RAW_REPCOMP_OUTDIR} && ${MARVEL_PATH}/bin/LAq${FIX_LAQ_OPT} -T trim0_d${RAW_FIX_LAQ_QTRIMCUTOFF}_s${RAW_FIX_LAQ_MINSEG}_repcomp -Q q0_d${RAW_FIX_LAQ_QTRIMCUTOFF}_s${RAW_FIX_LAQ_MINSEG}_repcomp ${RAW_DB%.db} -b ${x} ${RAW_DAZZ_DB%.db}.repcompFilt.${x}.las && cd ${myCWD}"
		done > fix_${sID}_LAq_block_${RAW_DB%.db}.${slurmID}.plan
    	echo "MARVEL LAq $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_${sID}_LAq_block_${RAW_DB%.db}.${slurmID}.version
	### 09_TKmerge    	                 
    elif [[ ${currentStep} -eq 9 ]]
    then
        ### clean up plans 
        for x in $(ls fix_${sID}_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done  
           
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
        echo "cd ${RAW_REPCOMP_OUTDIR} && ${MARVEL_PATH}/bin/TKmerge${FIX_TKMERGE_OPT} ${RAW_DB%.db} trim0_d${RAW_FIX_LAQ_QTRIMCUTOFF}_s${RAW_FIX_LAQ_MINSEG}_repcomp && cp .${RAW_DB%.db}.trim0_d${RAW_FIX_LAQ_QTRIMCUTOFF}_s${RAW_FIX_LAQ_MINSEG}_repcomp.a2 .${RAW_DB%.db}.trim0_d${RAW_FIX_LAQ_QTRIMCUTOFF}_s${RAW_FIX_LAQ_MINSEG}_repcomp.d2 ${myCWD}/ && cd ${myCWD}" > fix_${sID}_TKmerge_block_${RAW_DB%.db}.${slurmID}.plan
        echo "cd ${RAW_REPCOMP_OUTDIR} && ${MARVEL_PATH}/bin/TKmerge${FIX_TKMERGE_OPT} ${RAW_DB%.db} q0_d${RAW_FIX_LAQ_QTRIMCUTOFF}_s${RAW_FIX_LAQ_MINSEG}_repcomp&& cp .${RAW_DB%.db}.q0_d${RAW_FIX_LAQ_QTRIMCUTOFF}_s${RAW_FIX_LAQ_MINSEG}_repcomp.a2 .${RAW_DB%.db}.q0_d${RAW_FIX_LAQ_QTRIMCUTOFF}_s${RAW_FIX_LAQ_MINSEG}_repcomp.d2 ${myCWD}/ && cd ${myCWD}" >> fix_${sID}_TKmerge_block_${RAW_DB%.db}.${slurmID}.plan
        echo "MARVEL TKmerge $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_${sID}_TKmerge_block_${RAW_DB%.db}.${slurmID}.version               
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
elif [[ ${RAW_PATCH_TYPE} -eq 2 ]]
then 
	
	if [[ -z "${RAW_DACCORD_INDIR}" ]]
	then
		RAW_DACCORD_INDIR=${RAW_DALIGN_OUTDIR}	
	fi
	
	fsuffix="dalignFilt"
	if [[ "${RAW_DACCORD_INDIR}" == "${RAW_REPCOMP_OUTDIR}" ]]
	then
		fsuffix="repcompFilt"
	fi	
	
    ### create sub-directory and link relevant DB and Track files
    if [[ ${currentStep} -eq 1 ]]
    then
        ### clean up plans 
    for x in $(ls fix_${sID}_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done                 

    	echo "if [[ -d ${RAW_DACCORD_OUTDIR}_${RAW_DACCORD_INDIR} ]]; then mv ${RAW_DACCORD_OUTDIR}_${RAW_DACCORD_INDIR} ${RAW_DACCORD_OUTDIR}_${RAW_DACCORD_INDIR}_$(date '+%Y-%m-%d_%H-%M-%S'); fi && mkdir ${RAW_DACCORD_OUTDIR}_${RAW_DACCORD_INDIR} && ln -s -r .${RAW_DB%.db}.* ${RAW_DB%.db}.db .${RAW_DAZZ_DB%.db}.* ${RAW_DAZZ_DB%.db}.db ${RAW_DACCORD_OUTDIR}_${RAW_DACCORD_INDIR}" > fix_${sID}_createSubDir_single_${RAW_DB%.db}.${slurmID}.plan
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_${sID}_createSubDir_single_${RAW_DB%.db}.${slurmID}.version
 	### 02-lassort
	elif [[ ${currentStep} -eq 2 ]]
    then
        ### clean up plans 
    	for x in $(ls fix_${sID}_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 

		setlassortOptions
		
		for x in $(seq 1 ${nblocks})
        do
        	echo "${LASTOOLS_PATH}/bin/lassort ${FIX_LASSORT_OPT} ${RAW_DACCORD_OUTDIR}_${RAW_DACCORD_INDIR}/${RAW_DAZZ_DB%.db}.${fsuffix}Sort.${x}.las ${RAW_DACCORD_INDIR}/${RAW_DAZZ_DB%.db}.${fsuffix}.${x}.las"
		done > fix_${sID}_lassort_block_${RAW_DB%.db}.${slurmID}.plan    	         
        echo "LASTOOLS lassort $(git --git-dir=${DACCORD_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_${sID}_lassort_block_${RAW_DB%.db}.${slurmID}.version
    ### 03-computeIntrinsicQV
	elif [[ ${currentStep} -eq 3 ]]
    then
        ### clean up plans 
        for x in $(ls fix_${sID}_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 				
		
		for x in $(seq 1 ${nblocks})
        do
        	echo "cd ${RAW_DACCORD_OUTDIR}_${RAW_DACCORD_INDIR} && ln -s -f ${RAW_DAZZ_DB%.db}.${fsuffix}Sort.${x}.las ${RAW_DAZZ_DB%.db}.${x}.${fsuffix}Sort.las && ${DACCORD_PATH}/bin/computeintrinsicqv2 -d${RAW_COV} ${RAW_DAZZ_DB%.db}.db ${RAW_DAZZ_DB%.db}.${x}.${fsuffix}Sort.las && unlink ${RAW_DAZZ_DB%.db}.${x}.${fsuffix}Sort.las && cd ${myCWD}"
		done > fix_${sID}_computeintrinsicqv2_block_${RAW_DB%.db}.${slurmID}.plan    	         
        echo "DACCORD computeintrinsicqv2 $(git --git-dir=${DACCORD_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_${sID}_computeintrinsicqv2_block_${RAW_DB%.db}.${slurmID}.version
	### 04_Catrack
	elif [[ ${currentStep} -eq 4 ]]
    then
        ### clean up plans 
        for x in $(ls fix_${sID}_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        echo "cd ${RAW_DACCORD_OUTDIR}_${RAW_DACCORD_INDIR} && ${DAZZLER_PATH}/bin/Catrack -v -f -d ${RAW_DAZZ_DB%.db}.db inqual && cp .${RAW_DAZZ_DB%.db}.inqual.anno .${RAW_DAZZ_DB%.db}.inqual.data ${myCWD}/ && cd ${myCWD}" > fix_${sID}_Catrack_single_${RAW_DB%.db}.${slurmID}.plan
		echo "DAZZ_DB Catrack $(git --git-dir=${DAZZLER_SOURCE_PATH}/DAZZ_DB/.git rev-parse --short HEAD)" > fix_${sID}_Catrack_single_${RAW_DB%.db}.${slurmID}.version                
    ### 05_lasdetectsimplerepeats
    elif [[ ${currentStep} -eq 5 ]]
    then
        ### clean up plans 
        for x in $(ls fix_${sID}_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
                
        OPT=""
        if [[ -z "${RAW_FIX_LASDETECTSIMPLEREPEATS_ERATE}" ]]
        then 
        	RAW_FIX_LASDETECTSIMPLEREPEATS_ERATE=0.35
   	 	fi 
   	 	
   	 	OPT="${OPT} -d$((RAW_COV/2)) -e${RAW_FIX_LASDETECTSIMPLEREPEATS_ERATE}"
    	
        for x in $(seq 1 ${nblocks})
        do
        	echo "cd ${RAW_DACCORD_OUTDIR}_${RAW_DACCORD_INDIR} && ${DACCORD_PATH}/bin/lasdetectsimplerepeats ${OPT} ${RAW_DAZZ_DB%.db}.rep.${x}.data ${RAW_DAZZ_DB%.db}.db ${RAW_DAZZ_DB%.db}.${fsuffix}Sort.${x}.las && cd ${myCWD}"
		done > fix_${sID}_lasdetectsimplerepeats_block_${RAW_DB%.db}.${slurmID}.plan
      	echo "DACCORD lasdetectsimplerepeats $(git --git-dir=${DACCORD_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_${sID}_lasdetectsimplerepeats_block_${RAW_DB%.db}.${slurmID}.version
    ### 06_mergeAndSortRepeats
    elif [[ ${currentStep} -eq 6 ]]
    then
        ### clean up plans 
        for x in $(ls fix_${sID}_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        files="${RAW_DAZZ_DB%.db}.rep.[0-9].data"
		if [[ ${nblocks} -gt 9 ]]
		then
			files="${files} ${RAW_DAZZ_DB%.db}.rep.[0-9][0-9].data"
		fi
		if [[ ${nblocks} -gt 99 ]]
		then
			files="${files} ${RAW_DAZZ_DB%.db}.rep.[0-9][0-9][0-9].data"
		fi
		if [[ ${nblocks} -gt 999 ]]
		then
			files="${files} ${RAW_DAZZ_DB%.db}.rep.[0-9][0-9][0-9][0-9].data"
		fi
		if [[ ${nblocks} -gt 9999 ]]
		then
			files="${files} ${RAW_DAZZ_DB%.db}.rep.[0-9][0-9][0-9][0-9][0-9].data"
		fi
    	if [[ ${nblocks} -gt 99999 ]]
        then
    		(>&2 echo "fix_${sID}_mergeAndSortRepeats: more than 99999 db blocks are not supported!!!")
        	exit 1	
    	fi
    	## sanity check 
    	cd ${RAW_DACCORD_OUTDIR}_${RAW_DACCORD_INDIR} && if [[ $(ls ${files} | wc -l) -ne ${nblocks} ]]; then exit 1; fi && cd ${myCWD}
    	echo "cd ${RAW_DACCORD_OUTDIR}_${RAW_DACCORD_INDIR} && cat ${files} | ${DACCORD_PATH}/bin/repsort ${RAW_DAZZ_DB%.db}.db > ${RAW_DAZZ_DB%.db}.rep.data && cd ${myCWD}" >> fix_${sID}_mergeAndSortRepeats_single_${RAW_DB%.db}.${slurmID}.plan
        echo "DACCORD repsort $(git --git-dir=${DACCORD_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_${sID}_mergeAndSortRepeats_single_${RAW_DB%.db}.${slurmID}.version
    ### 07_lasfilteralignments 
    elif [[ ${currentStep} -eq 7 ]]
    then
        ### clean up plans 
        for x in $(ls fix_${sID}_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        OPT=""
        
        if [[ -z "${RAW_FIX_LASFILTERALIGNMENTS_ERATE}" ]]
        then 
        	RAW_FIX_LASFILTERALIGNMENTS_ERATE=0.35
   	 	fi 
   	 	
   	 	OPT="${OPT} -e${RAW_FIX_LASFILTERALIGNMENTS_ERATE}"
    	
        for x in $(seq 1 ${nblocks})
        do
        	echo "cd ${RAW_DACCORD_OUTDIR}_${RAW_DACCORD_INDIR} && ${DACCORD_PATH}/bin/lasfilteralignments ${OPT} ${RAW_DAZZ_DB%.db}.${fsuffix}SortFilt1.${x}.las ${RAW_DAZZ_DB%.db}.db ${RAW_DAZZ_DB%.db}.${fsuffix}Sort.${x}.las && cd ${myCWD}"
		done > fix_${sID}_lasfilteralignments_block_${RAW_DB%.db}.${slurmID}.plan
      	echo "DACCORD lasfilteralignments $(git --git-dir=${DACCORD_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_${sID}_lasfilteralignments_block_${RAW_DB%.db}.${slurmID}.version
    ### 08_mergesym2
    elif [[ ${currentStep} -eq 8 ]]
    then
        ### clean up plans 
        for x in $(ls fix_${sID}_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        echo "cd ${RAW_DACCORD_OUTDIR}_${RAW_DACCORD_INDIR} && ${DACCORD_PATH}/bin/mergesym2 ${RAW_DAZZ_DB%.db}.${fsuffix}SortFilt1.sym ${RAW_DAZZ_DB%.db}.db ${RAW_DAZZ_DB%.db}.${fsuffix}SortFilt1.*.las.sym && cd ${myCWD}" > fix_${sID}_mergesym2_single_${RAW_DB%.db}.${slurmID}.plan
        echo "cd ${RAW_DACCORD_OUTDIR}_${RAW_DACCORD_INDIR} && rm ${RAW_DAZZ_DB%.db}.${fsuffix}SortFilt1.*.las.sym && cd ${myCWD}" >> fix_${sID}_mergesym2_single_${RAW_DB%.db}.${slurmID}.plan
        echo "DACCORD mergesym2 $(git --git-dir=${DACCORD_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_${sID}_mergesym2_single_${RAW_DB%.db}.${slurmID}.version        
	### 09_filtersym
    elif [[ ${currentStep} -eq 9 ]]
    then
        ### clean up plans 
        for x in $(ls fix_${sID}_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        OPT=""        
        
		if [[ -z "${RAW_FILT_FILTERSYM_VERBOSE}" ]]
        then
        	RAW_FILT_FILTERSYM_VERBOSE=1
   	 	fi 
   	 	
   	 	if [[ -n "${RAW_FILT_FILTERSYM_VERBOSE}" && ${RAW_FILT_FILTERSYM_VERBOSE} != 0 ]]
        then
   	 		OPT="--verbose" 
   	 	fi
   	 	
   	 	for x in $(seq 1 ${nblocks})
        do
    		echo "cd ${RAW_DACCORD_OUTDIR}_${RAW_DACCORD_INDIR} && ${DACCORD_PATH}/bin/filtersym ${OPT} ${RAW_DAZZ_DB%.db}.${fsuffix}SortFilt1.${x}.las ${RAW_DAZZ_DB%.db}.${fsuffix}SortFilt1.sym" 
		done > fix_${sID}_filtsym_block_${RAW_DB%.db}.${slurmID}.plan
      	echo "DACCORD filtsym $(git --git-dir=${DACCORD_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_${sID}_filtsym_block_${RAW_DB%.db}.${slurmID}.version                 
   	### 10_lasfilteralignmentsborderrepeats
    elif [[ ${currentStep} -eq 10 ]]
    then
        ### clean up plans 
        for x in $(ls fix_10_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
                
		OPT=""
        
		if [[ -z "${RAW_FILT_LASFILTERALIGNMENTSBORDERREPEATS_THREADS}" ]]
        then
        	RAW_FILT_LASFILTERALIGNMENTSBORDERREPEATS_THREADS=8
   	 	fi 
   	 	
   	 	OPT="-t${RAW_FILT_LASFILTERALIGNMENTSBORDERREPEATS_THREADS}"
   	 	
   	 	if [[ -z "${RAW_FILT_LASFILTERALIGNMENTSBORDERREPEATS_ERATE}" ]]
        then
        	RAW_FILT_LASFILTERALIGNMENTSBORDERREPEATS_ERATE=0.35
   	 	fi 
   	 	
   	 	OPT="${OPT} -e${RAW_FILT_LASFILTERALIGNMENTSBORDERREPEATS_ERATE}"
   	 	for x in $(seq 1 ${nblocks})
        do
    		echo "cd ${RAW_DACCORD_OUTDIR}_${RAW_DACCORD_INDIR} && ${DACCORD_PATH}/bin/lasfilteralignmentsborderrepeats ${OPT} ${RAW_DAZZ_DB%.db}.${fsuffix}SortFilt2.${x}.las ${RAW_DAZZ_DB%.db}.db ${RAW_DAZZ_DB%.db}.rep.data ${RAW_DAZZ_DB%.db}.${fsuffix}SortFilt1.${x}.las && cd ${mxCWD}" 
		done > fix_10_lasfilteralignmentsborderrepeats_block_${RAW_DB%.db}.${slurmID}.plan
      	echo "DACCORD lasfilteralignmentsborderrepeats $(git --git-dir=${DACCORD_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_10_lasfilteralignmentsborderrepeats_block_${RAW_DB%.db}.${slurmID}.version
  	### 11_mergesym2
    elif [[ ${currentStep} -eq 11 ]]
    then
        ### clean up plans 
        for x in $(ls fix_${sID}_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        OPT=""        
        echo "cd ${RAW_DACCORD_OUTDIR}_${RAW_DACCORD_INDIR} && ${DACCORD_PATH}/bin/mergesym2 ${OPT} ${RAW_DAZZ_DB%.db}.${fsuffix}SortFilt2.sym ${RAW_DAZZ_DB%.db}.db ${RAW_DAZZ_DB%.db}.${fsuffix}SortFilt2.*.las.sym && cd ${myCWD}" > fix_${sID}_mergesym2_single_${RAW_DB%.db}.${slurmID}.plan
        echo "cd ${RAW_DACCORD_OUTDIR}_${RAW_DACCORD_INDIR} && rm ${RAW_DAZZ_DB%.db}.${fsuffix}SortFilt2.*.las.sym && cd ${myCWD}" >> fix_${sID}_mergesym2_single_${RAW_DB%.db}.${slurmID}.plan
        echo "DACCORD mergesym2 $(git --git-dir=${DACCORD_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_${sID}_mergesym2_single_${RAW_DB%.db}.${slurmID}.version        
	### 12_filtersym
    elif [[ ${currentStep} -eq 12 ]]
    then
        ### clean up plans 
        for x in $(ls fix_${sID}_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        OPT=""
        
		if [[ -z "${RAW_FILT_FILTERSYM_VERBOSE}" ]]
        then
        	RAW_FILT_FILTERSYM_VERBOSE=1
   	 	fi 
   	 	
   	 	if [[ -n "${RAW_FILT_FILTERSYM_VERBOSE}" && ${RAW_FILT_FILTERSYM_VERBOSE} != 0 ]]
        then
   	 		OPT="--verbose" 
   	 	fi
   	 	
   	 	for x in $(seq 1 ${nblocks})
        do
    		echo "cd ${RAW_DACCORD_OUTDIR}_${RAW_DACCORD_INDIR} && ${DACCORD_PATH}/bin/filtersym ${OPT} ${RAW_DAZZ_DB%.db}.${fsuffix}SortFilt2.${x}.las ${RAW_DAZZ_DB%.db}.${fsuffix}SortFilt2.sym && cd ${myCWD}" 
		done > fix_${sID}_filtsym_block_${RAW_DB%.db}.${slurmID}.plan
      	echo "DACCORD filtsym $(git --git-dir=${DACCORD_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_${sID}_filtsym_block_${RAW_DB%.db}.${slurmID}.version
    ### 13_filterchainsraw
    elif [[ ${currentStep} -eq 13 ]]
    then
        ### clean up plans 
        for x in $(ls fix_${sID}_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
                
        OPT=""
        setLAfilterOptions
		if [[ -z "${RAW_FILT_FILTERCHAINSRAW_LEN}" ]]
        then
        	RAW_FILT_FILTERCHAINSRAW_LEN=4000
   	 	fi 
   	 	
   	 	OPT="-l${RAW_FILT_FILTERCHAINSRAW_LEN}"
        for x in $(seq 1 ${nblocks})
        do
    		echo "cd ${RAW_DACCORD_OUTDIR}_${RAW_DACCORD_INDIR} && ${DACCORD_PATH}/bin/filterchainsraw ${OPT} ${RAW_DAZZ_DB%.db}.${fsuffix}SortFilt2Chain.${x}.las ${RAW_DAZZ_DB%.db}.db ${RAW_DAZZ_DB%.db}.${fsuffix}SortFilt2.${x}.las && ${MARVEL_PATH}/bin/LAfilter ${FIX_LAFILTER_OPT} ${RAW_DB%.db}.db ${RAW_DAZZ_DB%.db}.${fsuffix}SortFilt2Chain.${x}.las ${RAW_DAZZ_DB%.db}.${fsuffix}SortFilt2Chain2.${x}.las && cd ${myCWD}" 
		done > fix_${sID}_filterchainsraw_block_${RAW_DB%.db}.${slurmID}.plan
    	echo "DACCORD filterchainsraw $(git --git-dir=${DACCORD_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_${sID}_filterchainsraw_block_${RAW_DB%.db}.${slurmID}.version      	
    ### 14_daccord
    elif [[ ${currentStep} -eq 14 ]]
    then
        ### clean up plans 
    	for x in $(ls fix_${sID}_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
	
		setDaccordOptions        
		
		for x in $(seq 1 ${nblocks})
		do
    		echo "cd ${RAW_DACCORD_OUTDIR}_${RAW_DACCORD_INDIR} && ${DACCORD_PATH}/bin/daccord ${FIX_DACCORD_OPT} --eprofonly -E${RAW_DAZZ_DB%.db}.${fsuffix}SortFilt2Chain2.${x}.eprof ${RAW_DAZZ_DB%.db}.${fsuffix}SortFilt2Chain2.${x}.las ${RAW_DAZZ_DB%.db}.db && ${DACCORD_PATH}/bin/daccord ${FIX_DACCORD_OPT} -E${RAW_DAZZ_DB%.db}.${fsuffix}SortFilt2Chain2.${x}.eprof ${RAW_DAZZ_DB%.db}.${fsuffix}SortFilt2Chain2.${x}.las ${RAW_DAZZ_DB%.db}.db > ${RAW_DAZZ_DB%.db}.${fsuffix}SortFilt2Chain2.${x}.dac.fasta && cd ${myCWD}"
		done > fix_${sID}_daccord_block_${RAW_DB%.db}.${slurmID}.plan
        echo "DACCORD daccord $(git --git-dir=${DACCORD_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_${sID}_daccord_block_${RAW_DB%.db}.${slurmID}.version
   	### 15_computeextrinsicqv
    elif [[ ${currentStep} -eq 15 ]]
    then
        ### clean up plans 
    	for x in $(ls fix_${sID}_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        
        if [[ ${nblocks} -lt 10 ]]
		then
			files="${RAW_DAZZ_DB%.db}.${fsuffix}SortFilt2Chain2.[0-9].dac.fasta"
		elif [[ ${nblocks} -lt 100 ]]
		then
			files="${RAW_DAZZ_DB%.db}.${fsuffix}SortFilt2Chain2.[0-9].dac.fasta ${RAW_DAZZ_DB%.db}.${fsuffix}SortFilt2Chain2.[0-9][0-9].dac.fasta"
		elif [[ ${nblocks} -lt 1000 ]]
		then
			files="${RAW_DAZZ_DB%.db}.${fsuffix}SortFilt2Chain2.[0-9].dac.fasta ${RAW_DAZZ_DB%.db}.${fsuffix}SortFilt2Chain2.[0-9][0-9].dac.fasta ${RAW_DAZZ_DB%.db}.${fsuffix}SortFilt2Chain2.[0-9][0-9][0-9].dac.fasta"
		elif [[ ${nblocks} -lt 10000 ]]
		then
			files="${RAW_DAZZ_DB%.db}.${fsuffix}SortFilt2Chain2.[0-9].dac.fasta ${RAW_DAZZ_DB%.db}.${fsuffix}SortFilt2Chain2.[0-9][0-9].dac.fasta ${RAW_DAZZ_DB%.db}.${fsuffix}SortFilt2Chain2.[0-9][0-9][0-9].dac.fasta ${RAW_DAZZ_DB%.db}.${fsuffix}SortFilt2Chain2.[0-9][0-9][0-9][0-9].dac.fasta"
		else
    		(>&2 echo "fix_${sID}_computeextrinsicqv_single_${RAW_DB%.db}.${slurmID}.: more than 99999 db blocks are not supported!!!")
        	exit 1	
    	fi
    	
    	OPT=""
		if [[ -n "${RAW_FILT_COMPUTEEXTRINSICQ_THREADS}" ]]
        then
        	OPT="${OPT} -t${RAW_FILT_COMPUTEEXTRINSICQ_THREADS}"
   	 	fi
		echo "cd ${RAW_DACCORD_OUTDIR}_${RAW_DACCORD_INDIR} && cat ${files} > ${RAW_DAZZ_DB%.db}.${fsuffix}SortFilt2Chain2.dac.fasta && ${DACCORD_PATH}/bin/computeextrinsicqv${OPT} ${RAW_DAZZ_DB%.db}.${fsuffix}SortFilt2Chain2.dac.fasta ${RAW_DAZZ_DB%.db}.db && cd ${myCWD}" > fix_${sID}_computeextrinsicqv_single_${RAW_DB%.db}.${slurmID}.plan
        echo "DACCORD computeextrinsicqv $(git --git-dir=${DACCORD_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_${sID}_computeextrinsicqv_single_${RAW_DB%.db}.${slurmID}.version
    ### 16_split
    elif [[ ${currentStep} -eq 16 ]]
    then
        ### clean up plans 
        for x in $(ls fix_${sID}_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
		if [[ -z "${RAW_FIX_SPLIT_DIVIDEBLOCK}" ]]
		then 
			RAW_FIX_SPLIT_DIVIDEBLOCK=10
		fi
		
		# create folder structure
		for x in $(seq 1 ${nblocks})
		do
			directory="${RAW_DACCORD_OUTDIR}_${RAW_DACCORD_INDIR}/${RAW_FIX_SPLIT_TYPE}_s${x}"
			if [[ -d "${directory}" ]]
			then
					mv ${directory} ${directory}_$(date '+%Y-%m-%d_%H-%M-%S'); 
			fi
			
			mkdir -p ${directory}			
		done 
                
        setHaploSplitOptions
        
        for x in $(seq 1 ${nblocks})
		do
			for y in $(seq 0 $((RAW_FIX_SPLIT_DIVIDEBLOCK-1)))
			do
				echo "cd ${RAW_DACCORD_OUTDIR}_${RAW_DACCORD_INDIR} && ${DACCORD_PATH}/bin/${FIX_SPLIT_OPT} -E${RAW_DAZZ_DB%.db}.${fsuffix}SortFilt2Chain2.${x}.eprof -J${y},${RAW_FIX_SPLIT_DIVIDEBLOCK} ${RAW_FIX_SPLIT_TYPE}_s${x}/${RAW_DAZZ_DB%.db}.${fsuffix}SortFilt2Chain2Split.${y}.${x}.las ${RAW_DAZZ_DB%.db}.${fsuffix}SortFilt2Chain2.${x}.dac.fasta ${RAW_DAZZ_DB%.db}.${fsuffix}SortFilt2Chain2.${x}.las ${RAW_DAZZ_DB%.db}.db && cd ${myCWD}"		
			done	    		
		done > fix_${sID}_split_block_${RAW_DB%.db}.${slurmID}.plan
        echo "DACCORD ${RAW_FIX_SPLIT_TYPE} $(git --git-dir=${DACCORD_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_${sID}_split_block_${RAW_DB%.db}.${slurmID}.version
	### 17_LAmerge 
    elif [[ ${currentStep} -eq 17 ]]
    then
        ### clean up plans 
        for x in $(ls fix_${sID}_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        		        
        setLAmergeOptions
        setHaploSplitOptions
        
        for x in $(seq 1 ${nblocks})
		do
			echo "cd ${RAW_DACCORD_OUTDIR}_${RAW_DACCORD_INDIR} && ${MARVEL_PATH}/bin/LAmerge ${FIX_LAMERGE_OPT} ${RAW_DB%.db} ${RAW_DAZZ_DB%.db}.${fsuffix}SortFilt2Chain2_${RAW_FIX_SPLIT_TYPE}.${x}.keep.las ${RAW_FIX_SPLIT_TYPE}_s${x}/${RAW_DAZZ_DB%.db}.${fsuffix}SortFilt2Chain2Split.*.${x}.las ${myCWD}/identity/${RAW_DAZZ_DB%.db}.identity.${x}.las && cd ${myCWD}"
			echo "cd ${RAW_DACCORD_OUTDIR}_${RAW_DACCORD_INDIR} && ${MARVEL_PATH}/bin/LAmerge ${FIX_LAMERGE_OPT} ${RAW_DB%.db} ${RAW_DAZZ_DB%.db}.${fsuffix}SortFilt2Chain2_${RAW_FIX_SPLIT_TYPE}.${x}.drop.las ${RAW_FIX_SPLIT_TYPE}_s${x}/${RAW_DAZZ_DB%.db}.${fsuffix}SortFilt2Chain2Split.*.${x}_drop.las ${myCWD}/identity/${RAW_DAZZ_DB%.db}.identity.${x}.las && cd ${myCWD}"	
		done > fix_${sID}_LAmerge_block_${RAW_DB%.db}.${slurmID}.plan
        echo "MARVEL LAmerge $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_${sID}_LAmerge_block_${RAW_DB%.db}.${slurmID}.version
	### 18_LAfix    
    elif [[ ${currentStep} -eq 18 ]]
    then
        ### clean up plans 
        for x in $(ls fix_${sID}_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        ### find and set LAfix options
        
		if [[ "${fsuffix}" == "dalignFilt" ]]
		then
		   	setLAfixOptions dalign
		else
			setLAfixOptions repcomp
		fi
		
		setHaploSplitOptions
		
        mkdir -p ${RAW_FIX_LAFIX_PATH}_daccord_${RAW_FIX_SPLIT_TYPE}
		
		addopt=""

        for x in $(seq 1 ${nblocks})
        do 
        	if [[ -n ${RAW_FIX_LAFIX_TRIMFILEPREFIX} ]]
        	then 
        		addopt="-T${RAW_FIX_LAFIX_TRIMFILEPREFIX}_${x}.txt "
        	fi
            echo "${MARVEL_PATH}/bin/LAfix${FIX_LAFIX_OPT} ${addopt}${RAW_DB%.db} ${RAW_DACCORD_OUTDIR}_${RAW_DACCORD_INDIR}/${RAW_DAZZ_DB%.db}.${fsuffix}SortFilt2Chain2_${RAW_FIX_SPLIT_TYPE}.${x}.keep.las ${RAW_FIX_LAFIX_PATH}_daccord_${RAW_FIX_SPLIT_TYPE}/${RAW_DB%.db}.${x}${RAW_FIX_LAFIX_FILESUFFIX}.fasta"
    	done > fix_${sID}_LAfix_block_${RAW_DB%.db}.${slurmID}.plan
    echo "MARVEL LAfix $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_${sID}_LAfix_block_${RAW_DB%.db}.${slurmID}.version                
	else
        (>&2 echo "step ${currentStep} in FIX_FILT_TYPE ${FIX_FILT_TYPE} not supported")
        (>&2 echo "valid steps are: ${myTypes[${FIX_FILT_TYPE}]}")
        exit 1            
    fi
elif [[ ${RAW_PATCH_TYPE} -eq 3 ]]
then
  	if [[ ${currentStep} -eq 1 ]]
    then
        ### clean up plans 
        for x in $(ls fix_${sID}_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
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
	    	echo "${SUBMIT_SCRIPTS_PATH}/patchingStats.sh ${configFile} 1" > fix_${sID}_patchingStats_block_${RAW_DB%.db}.${slurmID}.plan
		fi
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > fix_${sID}_patchingStats_block_${RAW_DB%.db}.${slurmID}.version
    else
		(>&2 echo "step ${currentStep} in RAW_PATCH_TYPE ${RAW_PATCH_TYPE} not supported")
        (>&2 echo "valid steps are: ${myTypes[${RAW_PATCH_TYPE}]}")
        exit 1        
    fi
else
    (>&2echo "unknown RAW_PATCH_TYPE ${RAW_PATCH_TYPE}")    
    (>&2 echo "supported types")
    x=0; while [ $x -lt ${#myTypes[*]} ]; do (>&2 echo "${myTypes[${x}]}"); done
    exit 1
fi

exit 0
