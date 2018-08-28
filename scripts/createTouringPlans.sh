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

if [[ ! -n "${FIX_TOUR_TYPE}" ]]
then 
    (>&2 echo "cannot create touring scripts if variable FIX_TOUR_TYPE is not set.")
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

fixblocks=$(getNumOfDbBlocks ${FIX_DB%.db}.db)

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

function setOGbuildOptions()
{
    TOUR_OGBUILD_OPT=""
    
    if [[ -n ${FIX_TOUR_OGBUILD_CONT} && ${FIX_TOUR_OGBUILD_CONT} -ne 0 ]]
    then
        TOUR_OGBUILD_OPT="${TOUR_OGBUILD_OPT} -c ${FIX_TOUR_OGBUILD_CONT}"
    fi
    if [[ -n ${FIX_TOUR_OGBUILD_SPLIT} && ${FIX_TOUR_OGBUILD_SPLIT} -ne 0 ]]
    then
        TOUR_OGBUILD_OPT="${TOUR_OGBUILD_OPT} -s "
    fi
    if [[ -n ${FIX_TOUR_OGBUILD_TRIM} && ${FIX_TOUR_OGBUILD_TRIM} -ne 0 ]]
    then
        if [[ -z ${SCRUB_LAQ_OPT} ]]
        then 
            setLAqOptions
        fi

        if [ -n ${FIX_SCRUB_TYPE} ]
        then
            if [[ ${FIX_SCRUB_TYPE} -eq 1 ]]
            then 
                TOUR_OGBUILD_OPT="${TOUR_OGBUILD_OPT} -t trim1_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}_dalign"
            elif [[ ${FIX_SCRUB_TYPE} -eq 2 ]]
            then 
                TOUR_OGBUILD_OPT="${TOUR_OGBUILD_OPT} -t trim1_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}_repcomp"
            elif [[ ${FIX_SCRUB_TYPE} -eq 3 ]]
            then 
                TOUR_OGBUILD_OPT="${TOUR_OGBUILD_OPT} -t trim1_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}_forcealign"
            fi
        else
            TOUR_OGBUILD_OPT="${TOUR_OGBUILD_OPT} -t trim1_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}"
        fi
    fi
}

function setOGtourOptions()
{
    TOUR_OGTOUR_OPT=""
    if [[ -n ${FIX_TOUR_OGTOUR_CIRCULAR} && ${FIX_TOUR_OGTOUR_CIRCULAR} -ne 0 ]]
    then
        TOUR_OGTOUR_OPT="${TOUR_OGTOUR_OPT} -c"
    fi
    if [[ -n ${FIX_TOUR_OGTOUR_DROPINV} && ${FIX_TOUR_OGTOUR_DROPINV} -ne 0 ]]
    then
        TOUR_OGTOUR_OPT="${TOUR_OGTOUR_OPT} -d"
    fi    
    if [[ -n ${FIX_TOUR_OGTOUR_LOOKAHAED} && ${FIX_TOUR_OGTOUR_LOOKAHAED} -gt 0 ]]
    then
        TOUR_OGTOUR_OPT="${TOUR_OGTOUR_OPT} -l ${FIX_TOUR_OGTOUR_LOOKAHAED}"
    fi    
}


function setLAfilterOptions()
{
    FILT_LAFILTER_OPT=""

    if [[ -z ${FIX_FILT_OUTDIR} ]]
    then
        FIX_FILT_OUTDIR="filtered"
    fi

    if [[ -n ${FIX_SCRUB_TYPE} && ${FIX_SCRUB_TYPE} -eq 1 ]]
    then 
        if [[ $(echo ${FIX_FILT_OUTDIR} | awk -F _ '{print $NF}') != "dalign" ]]
        then
            FIX_FILT_OUTDIR="${FIX_FILT_OUTDIR}_dalign"
        fi
    elif [[ -n ${FIX_SCRUB_TYPE} && ${FIX_SCRUB_TYPE} -eq 2 ]]
    then 
        if [[ $(echo ${FIX_FILT_OUTDIR} | awk -F _ '{print $NF}') != "repcomp" ]]
        then
            FIX_FILT_OUTDIR="${FIX_FILT_OUTDIR}_repcomp"
        fi
    elif [[ -n ${FIX_SCRUB_TYPE} && ${FIX_SCRUB_TYPE} -eq 3 ]]
    then 
        if [[ $(echo ${FIX_FILT_OUTDIR} | awk -F _ '{print $NF}') != "forcealign" ]]
        then
            FIX_FILT_OUTDIR="${FIX_FILT_OUTDIR}_forcealign"
        fi
    fi
}

function settour2fastaOptions()
{
    TOUR_2FASTA_OPT=""
    if [[ -n ${FIX_TOUR_2FASTA_SPLIT} && ${FIX_TOUR_2FASTA_SPLIT} -ne 0 ]]
    then
        TOUR_2FASTA_OPT="${TOUR_2FASTA_OPT} -s"
    fi
    if [[ -n ${FIX_TOUR_2FASTA_TRIM} && ${FIX_TOUR_2FASTA_TRIM} -ne 0 ]]
    then
        if [[ -z ${SCRUB_LAQ_OPT} ]]
        then 
            setLAqOptions
        fi

        if [ -n ${FIX_SCRUB_TYPE} ]
        then
            if [[ ${FIX_SCRUB_TYPE} -eq 1 ]]
            then 
                TOUR_2FASTA_OPT="${TOUR_2FASTA_OPT} -t trim1_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}_dalign"
            elif [[ ${FIX_SCRUB_TYPE} -eq 2 ]]
            then 
                TOUR_2FASTA_OPT="${TOUR_2FASTA_OPT} -t trim1_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}_repcomp"
            elif [[ ${FIX_SCRUB_TYPE} -eq 3 ]]
            then 
                TOUR_2FASTA_OPT="${TOUR_2FASTA_OPT} -t trim1_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}_forcealign"
            fi
        else
            TOUR_2FASTA_OPT="${TOUR_2FASTA_OPT} -t trim1_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}"
        fi
    fi
}

function setOGlayoutOptions()
{
    TOUR_OGLAYOUT_OPT=""

    if [[ -n ${FIX_TOUR_OGLAYOUT_VERBOSE} && ${FIX_TOUR_OGLAYOUT_VERBOSE} -ne 0 ]]
    then
        TOUR_OGLAYOUT_OPT="${TOUR_OGLAYOUT_OPT} -v"
    fi
    if [[ -n ${FIX_TOUR_OGLAYOUT_DIST} && ${FIX_TOUR_OGLAYOUT_DIST} -ne 0 ]]
    then
        TOUR_OGLAYOUT_OPT="${TOUR_OGLAYOUT_OPT} -d ${FIX_TOUR_OGLAYOUT_DIST}"
    fi
    if [[ -n ${FIX_TOUR_OGLAYOUT_RMREVERSEEDGE} && ${FIX_TOUR_OGLAYOUT_RMREVERSEEDGE} -ne 0 ]]
    then
        TOUR_OGLAYOUT_OPT="${TOUR_OGLAYOUT_OPT} -R"
    fi
}

#type-0 steps: 1-OGbuild, 2-OGtour, 3-tour2fasta, 4-OGlayout, 5-marvelStats
if [[ ${FIX_TOUR_TYPE} -eq 0 ]]
then 
    ### OGbuild
    if [[ ${currentStep} -eq 1 ]]
    then
        ### clean up plans 
        for x in $(ls tour_01_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        
        setLAfilterOptions
        ### find and set OGbuild options 
        setOGbuildOptions
        ### create OGbuild commands
        echo "if [[ -d ${FIX_FILT_OUTDIR}/tour ]]; then rm -rf ${FIX_FILT_OUTDIR}/tour; fi" > tour_01_OGbuild_single_${FIX_DB%.db}.${slurmID}.plan
        echo "mkdir -p ${FIX_FILT_OUTDIR}/tour" >> tour_01_OGbuild_single_${FIX_DB%.db}.${slurmID}.plan
        echo "${MARVEL_PATH}/bin/OGbuild${TOUR_OGBUILD_OPT} ${FIX_FILT_OUTDIR}/${FIX_DB%.db} ${FIX_FILT_OUTDIR}/${FIX_DB%.db}.filt.las ${FIX_FILT_OUTDIR}/tour/${FIX_DB%.db}.filt" >> tour_01_OGbuild_single_${FIX_DB%.db}.${slurmID}.plan
    ### OGtour
    elif [[ ${currentStep} -eq 2 ]]
    then
        ### clean up plans 
        for x in $(ls tour_02_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        
        if [[ -z ${FILT_LAFILTER_OPT} ]]
        then
            setLAfilterOptions
        fi
        ### find and set OGbuild options 
        setOGtourOptions
        ### create OGbuild commands    
        for x in ${FIX_FILT_OUTDIR}/tour/*[0-9].graphml; 
        do 
            if [[ -s ${x} ]]
            then
                echo "${MARVEL_PATH}/scripts/OGtour.py${TOUR_OGTOUR_OPT} ${FIX_FILT_OUTDIR}/${FIX_DB} $x"
            fi 
        done > tour_02_OGtour_block_${FIX_DB%.db}.${slurmID}.plan        
    ### tour2fasta
    elif [[ ${currentStep} -eq 3 ]]
    then
        ### clean up plans 
        for x in $(ls tour_03_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        
        if [[ -z ${FILT_LAFILTER_OPT} ]]
        then
            setLAfilterOptions
        fi
        ### find and set OGbuild options 
        settour2fastaOptions
        
        for x in ${FIX_FILT_OUTDIR}/tour/*[0-9].tour.paths;
        do 
            if [[ -s ${x} ]]
            then
                echo "${MARVEL_PATH}/scripts/tour2fasta.py${TOUR_2FASTA_OPT} ${FIX_FILT_OUTDIR}/${FIX_DB} ${x%.tour.paths}.graphml $x"
            fi
        done > tour_03_tour2fasta_block_${FIX_DB%.db}.${slurmID}.plan
    ### OGlayout
    elif [[ ${currentStep} -eq 4 ]]
    then
        ### clean up plans 
        for x in $(ls tour_04_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        
        if [[ -z ${FILT_LAFILTER_OPT} ]]
        then
            setLAfilterOptions
        fi
        ### find and set OGbuild options 
        setOGlayoutOptions   

        for x in ${FIX_FILT_OUTDIR}/tour/*[0-9].tour.paths; 
        do 
            if [[ -s ${x} ]]
            then
                echo "${MARVEL_PATH}/bin/OGlayout${TOUR_OGLAYOUT_OPT} ${x%.paths}.graphml ${x%.paths}.layout.dot" 
            fi
        done > tour_04_OGlayout_block_${FIX_DB%.db}.${slurmID}.plan
    ### marvelStats
    elif [[ ${currentStep} -eq 5 ]]
    then
        ### clean up plans 
        for x in $(ls tour_05_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        
        if [[ -z ${FILT_LAFILTER_OPT} ]]
        then
            setLAfilterOptions
        fi
        ### run slurm stats - on the master node !!! Because sacct is not available on compute nodes
        bash ${SUBMIT_SCRIPTS_PATH}/slurmStats.sh ${configFile} > ${FIX_FILT_OUTDIR}/stats.${FIX_FILT_OUTDIR}.log
        ### create assemblyStats plan 
        echo "${SUBMIT_SCRIPTS_PATH}/assemblyStats.sh ${configFile} >> ${FIX_FILT_OUTDIR}/stats.${FIX_FILT_OUTDIR}.log" > tour_05_marvelStats_block_${FIX_DB%.db}.${slurmID}.plan
    else
        (>&2 echo "step ${currentStep} in FIX_TOUR_TYPE ${FIX_TOUR_TYPE} not supported")
        (>&2 echo "valid steps are: #type-0 steps: 1-OGbuild, 2-OGtour, 3-tour2fasta, 4-OGlayout, 5-marvelStats")
        exit 1            
    fi
else
    echo "unknown FIX_TOUR_TYPE ${FIX_TOUR_TYPE}"
    exit 1
fi

exit 0