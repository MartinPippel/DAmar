#!/bin/bash

config=$1
faFolder=$2

if [[ ! -f ${config} ]]
then 
  echo "config ${config} not available"
  exit 1
fi

if [[ ! -d ${faFolder} ]]
then 
  echo "dir ${faFolder} not available"
  exit 1
fi

source ${config}



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
    
    for x in $(seq 0 $((${numRepeatTracks}-1)))
    do 
        tmp=""
        tmp="${tmp} -l ${FIX_SCRUB_LAREPEAT_LEAVE_COV[$x]}"
        tmp="${tmp} -h ${FIX_SCRUB_LAREPEAT_ENTER_COV[$x]}"
        if [[ ${FIX_SCRUB_LAREPEAT_COV[$x]} -ne -1 ]]
        then 
            tmp="${tmp} -c ${FIX_SCRUB_LAREPEAT_COV[$x]}"
            tmp="${tmp} -t repeats_c${FIX_SCRUB_LAREPEAT_COV[$x]}_l${FIX_SCRUB_LAREPEAT_LEAVE_COV[$x]}h${FIX_SCRUB_LAREPEAT_ENTER_COV[$x]}"
        else
            tmp="${tmp} -t repeats_calCov_l${FIX_SCRUB_LAREPEAT_LEAVE_COV[$x]}h${FIX_SCRUB_LAREPEAT_ENTER_COV[$x]}"
        fi
        SCRUB_LAREPEAT_OPT[$x]=${tmp}
    done 
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

function setLAfilterOptions()
{
    FILT_LAFILTER_OPT=""

    if [[ -z ${FIX_FILT_FSUFFIX} ]]
    then
        FIX_FILT_FSUFFIX=filtered
    fi

    if [[ -n ${FIX_FILT_LAFILTER_NREP} && ${FIX_FILT_LAFILTER_NREP} -ne 0 ]]
    then
        FILT_LAFILTER_OPT="${FILT_LAFILTER_OPT} -n ${FIX_FILT_LAFILTER_NREP}"
    fi
    if [[ -n ${FIX_FILT_LAFILTER_VERBOSE} && ${FIX_FILT_LAFILTER_VERBOSE} -ne 0 ]]
    then
        FILT_LAFILTER_OPT="${FILT_LAFILTER_OPT} -v"
    fi
    if [[ -n ${FIX_FILT_LAFILTER_PURGE} && ${FIX_FILT_LAFILTER_PURGE} -ne 0 ]]
    then
        FILT_LAFILTER_OPT="${FILT_LAFILTER_OPT} -p"
    fi
    if [[ -n ${FIX_FILT_LAFILTER_OLEN} && ${FIX_FILT_LAFILTER_OLEN} -ne 0 ]]
    then
        FILT_LAFILTER_OPT="${FILT_LAFILTER_OPT} -o ${FIX_FILT_LAFILTER_OLEN}"
    fi    
    if [[ -n ${FIX_FILT_LAFILTER_RLEN} && ${FIX_FILT_LAFILTER_RLEN} -ne 0 ]]
    then
        FILT_LAFILTER_OPT="${FILT_LAFILTER_OPT} -l ${FIX_FILT_LAFILTER_RLEN}"
    fi   

    if [[ -n ${FIX_FILT_LAFILTER_REPEAT_IDX} ]]
    then 
        if [[ -z ${SCRUB_LAREPEAT_OPT} ]]
        then 
            setLArepeatOptions
        fi

        if [[ ${numRepeatTracks} -eq 0 || $((${FIX_FILT_LAFILTER_REPEAT_IDX}+1)) -gt ${#SCRUB_LAREPEAT_OPT[*]} ]]
        then 
            exit 1
        fi

        tmp=$(echo ${SCRUB_LAREPEAT_OPT[${FIX_FILT_LAFILTER_REPEAT_IDX}]} | awk '{print $NF}')_${RAW_REPMASK_LAREPEAT_REPEATTRACK}
        FIX_FILT_LAFILTER_REPEATTRACK=f${tmp}_${FIX_SCRUB_TANMASK_TRACK}_dust

        FILT_LAFILTER_OPT="${FILT_LAFILTER_OPT} -r ${FIX_FILT_LAFILTER_REPEATTRACK}"
    fi

    if [[ -n ${FIX_FILT_LAFILTER_DIF} && ${FIX_FILT_LAFILTER_DIF} -ne 0 ]]
    then
        FILT_LAFILTER_OPT="${FILT_LAFILTER_OPT} -d ${FIX_FILT_LAFILTER_DIF}"
    fi
    if [[ -n ${FIX_FILT_LAFILTER_TRIM} && ${FIX_FILT_LAFILTER_TRIM} -ne 0 ]]
    then
        if [[ -z ${SCRUB_LAQ_OPT} ]]
        then 
            setLAqOptions
        fi
        
        FILT_LAFILTER_OPT="${FILT_LAFILTER_OPT} -t trim1_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG} -T"
    fi
    if [[ -n ${FIX_FILT_LAFILTER_UBAS} ]]
    then
        FILT_LAFILTER_OPT="${FILT_LAFILTER_OPT} -u ${FIX_FILT_LAFILTER_UBAS}"
    fi
    if [[ -n ${FIX_FILT_LAFILTER_PRELOAD} && ${FIX_FILT_LAFILTER_PRELOAD} -ne 0 ]]
    then
        FILT_LAFILTER_OPT="${FILT_LAFILTER_OPT} -L"
    fi    
    if [[ -n ${FIX_FILT_LAFILTER_MERGEREPEATS} && ${FIX_FILT_LAFILTER_MERGEREPEATS} -ne 0 ]]
    then
        FILT_LAFILTER_OPT="${FILT_LAFILTER_OPT} -y ${FIX_FILT_LAFILTER_MERGEREPEATS}"
    fi    
    if [[ -n ${FIX_FILT_LAFILTER_MERGEREPEATTIPS} && ${FIX_FILT_LAFILTER_MERGEREPEATTIPS} -ne 0 ]]
    then
        FILT_LAFILTER_OPT="${FILT_LAFILTER_OPT} -Y ${FIX_FILT_LAFILTER_MERGEREPEATTIPS}"
    fi    
    if [[ -n ${FIX_FILT_LAFILTER_MINTIPCOV} && ${FIX_FILT_LAFILTER_MINTIPCOV} -ge 0 ]]
    then
        FILT_LAFILTER_OPT="${FILT_LAFILTER_OPT} -z ${FIX_FILT_LAFILTER_MINTIPCOV}"
    fi            
    if [[ -n ${FIX_FILT_LAFILTER_EXCLUDEREADS} ]]
    then
        FILT_LAFILTER_OPT="${FILT_LAFILTER_OPT} -x ${FIX_FILT_LAFILTER_EXCLUDEREADS}"
    fi    
    if [[ -n ${FIX_FILT_LAFILTER_RESOLVE_REPEATS} && ${FIX_FILT_LAFILTER_RESOLVE_REPEATS} -gt 0 && ${FIX_FILT_LAFILTER_RESOLVE_REPEATS} -lt 4 ]]
    then
        tmp=""
        mode="m"
        if [[ -n ${FIX_FILT_LAFILTER_RESOLVE_REPEATS_AGG} && ${FIX_FILT_LAFILTER_RESOLVE_REPEATS_AGG} -ne 0 ]]
        then
            mode="M"
        fi
        for x in $(seq 1 ${FIX_FILT_LAFILTER_RESOLVE_REPEATS})
        do
            tmp="${tmp}${mode}"
        done

        if [[ -z ${FIX_COV} ]]
        then 
            (>&2 echo "If FIX_FILT_LAFILTER_RESOLVE_REPEATS is set, then the FIX_COV variable has to be set with the apropriate coverage of the patched database ${FIX_DB%.db}.db!")
            exit 1
        fi 
        FILT_LAFILTER_OPT="${FILT_LAFILTER_OPT} -${tmp} ${FIX_COV}"
    fi    
}


setLAfilterOptions
setLAqOptions

frep=${FIX_FILT_LAFILTER_REPEATTRACK}
trim0=trim0_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}
trim1=trim1_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}

db=${DB%.db}
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

echo "#RUN OPTIONS"
echo "QUAL: ${SCRUB_LAQ1_QUAL}"
echo "FILT_REPEATS: ${FIX_SCRUB_LAQ_QTRIMCUTOFF}"
if [[ -n ${FIX_SCRUB_LAGAP_STITCH} ]]
then 
	echo "STITCH: ${FIX_SCRUB_LAGAP_STITCH}"
else
	echo "STITCH: NOT SET"
fi

if [[ -n ${FIX_FILT_LAFILTER_DIF} ]]
then
	echo "FIX_FILT_LAFILTER_DIF ${FIX_FILT_LAFILTER_DIF}"
else
	echo "FIX_FILT_LAFILTER_DIF NOT SET"
fi

if [[ -n ${FIX_FILT_LAFILTER_NREP} ]]
then
        echo "FIX_FILT_LAFILTER_NREP ${FIX_FILT_LAFILTER_NREP}"
else
        echo "FIX_FILT_LAFILTER_NREP NOT SET"
fi

if [[ -n ${FIX_FILT_LAFILTER_OLEN} ]]
then
        echo "FIX_FILT_LAFILTER_OLEN ${FIX_FILT_LAFILTER_OLEN}"
else
        echo "FIX_FILT_LAFILTER_OLEN NOT SET"
fi

if [[ -n ${FIX_FILT_LAFILTER_RLEN} ]]
then
        echo "FIX_FILT_LAFILTER_RLEN ${FIX_FILT_LAFILTER_RLEN}"
else
        echo "FIX_FILT_LAFILTER_RLEN NOT SET"
fi

if [[ -n ${FIX_TOUR_OGTOUR_CIRCULAR} ]]
then
        echo "FIX_TOUR_OGTOUR_CIRCULAR ${FIX_TOUR_OGTOUR_CIRCULAR}"
else
        echo "FIX_TOUR_OGTOUR_CIRCULAR NOT SET"
fi

if [[ -n ${FIX_TOUR_OGTOUR_LOOKAHAED} ]]
then
        echo "FIX_TOUR_OGTOUR_LOOKAHAED ${FIX_TOUR_OGTOUR_LOOKAHAED}"
else
        echo "FIX_TOUR_OGTOUR_LOOKAHAED NOT SET"
fi

if [[ -n ${FIX_TOUR_OGTOUR_DROPINV} ]]
then
        echo "FIX_TOUR_OGTOUR_DROPINV ${FIX_TOUR_OGTOUR_DROPINV}"
else
        echo "FIX_TOUR_OGTOUR_DROPINV NOT SET"
fi

cat ${faFolder}/*.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize}
${MARVEL_PATH}/bin/TKshow -s 4 -S -b 2000 ${db} $frep 
${MARVEL_PATH}/bin/TKshow -s 4 -S -b 2000 ${db} ${trim0}
${MARVEL_PATH}/bin/TKshow -s 4 -S -b 2000 ${db} ${trim1} 
cat $config
