#!/bin/bash 

configFile=$1
currentStep=$2
slurmID=$3

cwd=$(pwd)
echo "createRepmaskPlans2.sh config: ${configFile} currentStep: ${currentStep} ID: ${slurmID}"
echo "createRepmaskPlans2.sh cwd ${cwd}" 

if [[ ! -f ${configFile} ]]
then 
    (>&2 echo "cannot access config file ${configFile}")
    exit 1
fi

source ${configFile}

if [[ ! -n "${FIX_REPMASK_TYPE}" ]]
then 
    (>&2 echo "cannot create repmask jobs if varibale FIX_REPMASK_TYPE is not set.")
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
fixblocks=0

if [[ -f ${RAW_DB%.db}.db ]]
then 
	rawblocks=$(getNumOfDbBlocks ${RAW_DB%.db}.db)	
fi

if [[ -f ${FIX_DB%.db}.db ]]
then 
	fixblocks=$(getNumOfDbBlocks ${FIX_DB%.db}.db)	
fi

function setDBdustOptions()
{
    REPMASK_DBDUST_OPT=""
    if [[ -n ${FIX_REPMASK_DBDUST_BIAS} && ${FIX_REPMASK_DBDUST_BIAS} -ge 1 ]]
    then
        REPMASK_DBDUST_OPT="${REPMASK_DBDUST_OPT} -b"
    fi
}

function setCatrackOptions()
{
    REPMASK_CATRACK_OPT=""
    if [[ -n ${FIX_REPMASK_CATRACK_VERBOSE} && ${FIX_REPMASK_CATRACK_VERBOSE} -ge 1 ]]
    then
        REPMASK_CATRACK_OPT="${REPMASK_CATRACK_OPT} -v"
    fi
    if [[ -n ${FIX_REPMASK_CATRACK_DELETE} && ${FIX_REPMASK_CATRACK_DELETE} -ge 1 ]]
    then
        REPMASK_CATRACK_OPT="${REPMASK_CATRACK_OPT} -d"
    fi
    if [[ -n ${FIX_REPMASK_CATRACK_OVERWRITE} && ${FIX_REPMASK_CATRACK_OVERWRITE} -ge 1 ]]
    then
        REPMASK_CATRACK_OPT="${REPMASK_CATRACK_OPT} -f"
    fi
}

function setTANmaskOptions()
{
    REPMASK_TANMASK_OPT=""
    if [[ -n ${FIX_REPMASK_TANMASK_VERBOSE} && ${FIX_REPMASK_TANMASK_VERBOSE} -ge 1 ]]
    then
        REPMASK_TANMASK_OPT="${REPMASK_TANMASK_OPT} -v"
    fi
    if [[ -n ${FIX_REPMASK_TANMASK_MINLEN} && ${FIX_REPMASK_TANMASK_MINLEN} -ge 1 ]]
    then
        REPMASK_TANMASK_OPT="${REPMASK_TANMASK_OPT} -l${FIX_REPMASK_TANMASK_MINLEN}"
    fi
    if [[ -n ${FIX_REPMASK_TANMASK_TRACK} ]]
    then
        REPMASK_TANMASK_OPT="${REPMASK_TANMASK_OPT} -n${FIX_REPMASK_TANMASK_TRACK}"
    fi
}

function setDalignerOptions()
{
    REPMASK_DALIGNER_OPT=""
    if [[ -n ${FIX_REPMASK_DALIGNER_IDENTITY_OVLS} && ${FIX_REPMASK_DALIGNER_IDENTITY_OVLS} -gt 0 ]]
    then
        REPMASK_DALIGNER_OPT="${REPMASK_DALIGNER_OPT} -I"
    fi
    if [[ -n ${FIX_REPMASK_DALIGNER_KMER} && ${FIX_REPMASK_DALIGNER_KMER} -gt 0 ]]
    then
        REPMASK_DALIGNER_OPT="${REPMASK_DALIGNER_OPT} -k${FIX_REPMASK_DALIGNER_KMER}"
    fi
    if [[ -n ${FIX_REPMASK_DALIGNER_ERR} ]]
    then
        REPMASK_DALIGNER_OPT="${REPMASK_DALIGNER_OPT} -e${FIX_REPMASK_DALIGNER_ERR}"
    fi
    if [[ -n ${FIX_REPMASK_DALIGNER_OLEN} ]]
    then
        REPMASK_DALIGNER_OPT="${REPMASK_DALIGNER_OPT} -l${FIX_REPMASK_DALIGNER_OLEN}"
    fi    
    if [[ -n ${FIX_REPMASK_DALIGNER_MEM} && ${FIX_REPMASK_DALIGNER_MEM} -gt 0 ]]
    then
        REPMASK_DALIGNER_OPT="${REPMASK_DALIGNER_OPT} -M${FIX_REPMASK_DALIGNER_MEM}"
    fi    
    if [[ -n ${FIX_REPMASK_DALIGNER_HITS} ]]
    then
        REPMASK_DALIGNER_OPT="${REPMASK_DALIGNER_OPT} -h${FIX_REPMASK_DALIGNER_HITS}"
    fi        
    if [[ -n ${FIX_REPMASK_DALIGNER_T} ]]
    then
        REPMASK_DALIGNER_OPT="${REPMASK_DALIGNER_OPT} -t${FIX_REPMASK_DALIGNER_T}"
    fi  
    if [[ -n ${FIX_REPMASK_DALIGNER_MASK} ]]
    then
        for x in ${FIX_REPMASK_DALIGNER_MASK}
        do 
            REPMASK_DALIGNER_OPT="${REPMASK_DALIGNER_OPT} -m${x}"
        done
    fi
    if [[ -n ${FIX_REPMASK_DALIGNER_TRACESPACE} && ${FIX_REPMASK_DALIGNER_TRACESPACE} -gt 0 ]]
    then
        REPMASK_DALIGNER_OPT="${REPMASK_DALIGNER_OPT} -s${FIX_REPMASK_DALIGNER_TRACESPACE}"
    fi
    if [[ -n ${THREADS_daligner} ]]
    then 
        REPMASK_DALIGNER_OPT="${REPMASK_DALIGNER_OPT} -T${THREADS_daligner}"
    fi
}

function setLArepeatOptions()
{
    idx=$1
    REPMASK_LAREPEAT_OPT=""
    if [[ -n ${FIX_REPMASK_LAREPEAT_LOW} ]]
    then
        REPMASK_LAREPEAT_OPT="${REPMASK_LAREPEAT_OPT} -l ${FIX_REPMASK_LAREPEAT_LOW}"
    fi
    if [[ -n ${FIX_REPMASK_LAREPEAT_HGH} ]]
    then
        REPMASK_LAREPEAT_OPT="${REPMASK_LAREPEAT_OPT} -h ${FIX_REPMASK_LAREPEAT_HGH}"
    fi
    if [[ -n ${FIX_REPMASK_LAREPEAT_OLEN} ]]
    then
        REPMASK_LAREPEAT_OPT="${REPMASK_LAREPEAT_OPT} -o ${FIX_REPMASK_LAREPEAT_OLEN}"
    fi
    if [[ -n ${FIX_REPMASK_LAREPEAT_REPEATTRACK} ]]
    then
        if [[ -z ${idx} ]]
        then
          FIX_REPMASK_REPEATTRACK=${FIX_REPMASK_LAREPEAT_REPEATTRACK}
          REPMASK_LAREPEAT_OPT="${REPMASK_LAREPEAT_OPT} -t ${FIX_REPMASK_REPEATTRACK}"
        else 
            if [[ ${#FIX_REPMASK_LAREPEAT_COV[*]} -lt ${idx} ]]
            then 
                (>&2 echo "FIX_REPMASK_LAREPEAT_COV has lower the ${idx} elements")
                exit 1
            elif [[ ${#FIX_REPMASK_BLOCKCMP[*]} -lt ${idx} ]]
            then 
                (>&2 echo "FIX_REPMASK_BLOCKCMP has lower the ${idx} elements")
                exit 1
            fi
            FIX_REPMASK_REPEATTRACK=${FIX_REPMASK_LAREPEAT_REPEATTRACK}_B${FIX_REPMASK_BLOCKCMP[${idx}]}C${FIX_REPMASK_LAREPEAT_COV[${idx}]}
            REPMASK_LAREPEAT_OPT="${REPMASK_LAREPEAT_OPT} -t ${FIX_REPMASK_REPEATTRACK}"
        fi
    fi
    if [[ -n ${FIX_REPMASK_LAREPEAT_COV} ]]
    then
        if [[ -z ${idx} ]]
        then
            REPMASK_LAREPEAT_OPT="${REPMASK_LAREPEAT_OPT} -c ${FIX_REPMASK_LAREPEAT_COV}"
        else
            if [[ ${#FIX_REPMASK_LAREPEAT_COV[*]} -lt ${idx} ]]
            then 
                (>&2 echo "FIX_REPMASK_LAREPEAT_COV has lower the ${idx} elements")
                exit 1
            fi
            REPMASK_LAREPEAT_OPT="${REPMASK_LAREPEAT_OPT} -c ${FIX_REPMASK_LAREPEAT_COV[${idx}]}"
        fi 
    fi
}

function setTKmergeOptions()
{
    REPMASK_TKMERGE_OPT=""
    if [[ -n ${FIX_REPMASK_TKMERGE_DELETE} && ${FIX_REPMASK_TKMERGE_DELETE} -ge 1 ]]
    then
        REPMASK_TKMERGE_OPT="${REPMASK_TKMERGE_OPT} -d"
    fi
    if [ ! -n ${FIX_REPMASK_LAREPEAT_REPEATTRACK} ] ### fall back to default value!!!
    then
        FIX_REPMASK_LAREPEAT_REPEATTRACK="repeats"
    fi
}

function setDatanderOptions()
{
    ### find and set datander options 
    REPMASK_DATANDER_OPT=""
    if [[ -n ${FIX_REPMASK_DATANDER_THREADS} ]]
    then
        REPMASK_DATANDER_OPT="${REPMASK_DATANDER_OPT} -T${FIX_REPMASK_DATANDER_THREADS}"
    fi
    if [[ -n ${FIX_REPMASK_DATANDER_MINLEN} ]]
    then
        REPMASK_DATANDER_OPT="${REPMASK_DATANDER_OPT} -l${FIX_REPMASK_DATANDER_MINLEN}"
    fi
}

function setDBsplitOptions()
{
    SCRUB_DBSPLIT_OPT=""

    FIX_REPMASK_DBSPLIT_S=$(grep size ${RAW_DB%.db}.db | awk '{print $3}')

    if [[ -z ${FIX_REPMASK_DBSPLIT_S} ]]
    then
        (>&2 echo "database $db has not been partitioned. Run DBsplit first!")
        exit 1
    fi

    SCRUB_DBSPLIT_OPT="${SCRUB_DBSPLIT_OPT} -s${FIX_REPMASK_DBSPLIT_S}"
}


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

if [[ -z "${LASTOOLS_PATH}" ]]
then 
    (>&2 echo "ERROR - You have to set LASTOOLS_PATH. ")
    exit 1
fi

if [[ -z "${FIX_REPMASK_OUTDIR}" ]]
then
	FIX_REPMASK_OUTDIR=repmask	
fi

if [[ ${currentStep} -lt 10 ]]
then 
	sID=0${currentStep}
else
	sID=${currentStep}
fi
myCWD=$(pwd)

myTypes=("1-createFIX_DB, 2-DBdust, 3-Catrack, 4-datander, 5-TANmask, 6-Catrack, 7-daligner, 8-LAmerge, 9-LArepeat, 10-TKmerge, 11-daligner, 12-LAmerge, 13-LArepeat, 14-TKmerge")
# type_0 - steps: 1-createFIX_DB, 2-DBdust, 3-Catrack, 4-datander, 5-TANmask, 6-Catrack, 7-daligner, 8-LAmerge, 9-LArepeat, 10-TKmerge, 11-daligner, 12-LAmerge, 13-LArepeat, 14-TKmerge
if [[ ${FIX_REPMASK_TYPE} -eq 0 ]]
then	
	if [[ ${currentStep} -eq 1 ]]
    then
    	### clean up plans 
        for x in $(ls mask_${sID}_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        if [[ -f ${FIX_DB%.db}.db ]]; 
        then 
            (>&2 echo "p3_s1 DB ${FIX_DB%.db}.db already exists !!!")
            exit 1            
        fi
        
        if [[ -f ${FIX_DAZZ_DB%.db}.db ]]; 
        then 
            (>&2 echo "p3_s1 DAZZ DB ${FIX_DAZZ_DB%.db}.db already exists !!!")
            exit 1                        
        fi

        # get directory of patched reads
        setDBsplitOptions

        if [[ ! -f ${FIX_REPMASK_USELAFIX_PATH}/${RAW_DB%.db}.1${RAW_FIX_LAFIX_FILESUFFIX}.fasta ]]
        then
        (>&2 echo "cannot find file ${FIX_REPMASK_USELAFIX_PATH}/${RAW_DB%.db}.1${RAW_FIX_LAFIX_FILESUFFIX}.fasta")
            exit 1 
        fi 

        FIX_REPMASK_REPEATTRACK=""
        for x in $(seq 1 ${#FIX_REPMASK_BLOCKCMP[*]})
        do
            idx=$(($x-1))
            FIX_REPMASK_REPEATTRACK="${FIX_REPMASK_REPEATTRACK} -c ${FIX_REPMASK_LAREPEAT_REPEATTRACK}_B${FIX_REPMASK_BLOCKCMP[${idx}]}C${FIX_REPMASK_LAREPEAT_COV[${idx}]}"
        done
        
        for x in $(seq 1 ${rawblocks})
        do  
            if [[ ! -f ${FIX_REPMASK_USELAFIX_PATH}/${RAW_DB%.db}.${x}${RAW_FIX_LAFIX_FILESUFFIX}.fasta ]]
            then 
                (>&2 echo "cannot find file ${FIX_REPMASK_USELAFIX_PATH}/${RAW_DB%.db}.${x}${RAW_FIX_LAFIX_FILESUFFIX}.fasta")
                exit 1
            fi
            
            if [[ $x -eq 1 ]]
            then
            	echo "${MARVEL_PATH}/bin/FA2db -c source ${FIX_REPMASK_REPEATTRACK} -v ${FIX_DB%.db}.db ${FIX_REPMASK_USELAFIX_PATH}/${RAW_DB%.db}.1${RAW_FIX_LAFIX_FILESUFFIX}.fasta"
        		echo "${MARVEL_PATH}/bin/DBsplit${SCRUB_DBSPLIT_OPT} ${FIX_DB%.db}.db"                    	 
         	else
            	echo "${MARVEL_PATH}/bin/FA2db -v -c source ${FIX_REPMASK_REPEATTRACK} ${FIX_DB%.db}.db ${FIX_REPMASK_USELAFIX_PATH}/${RAW_DB%.db}.${x}${RAW_FIX_LAFIX_FILESUFFIX}.fasta"
        	fi
        done > mask_${sID}_createDB_single_${FIX_DB%.db}.${slurmID}.plan        

        ### convert corrected reads into a proper dazzler read format and create a valid dazzler db, that where read IDS do exactly map to the corressponding marvel db

        if [[ -d ${RAW_FIX_LAFIX_PATH}${ptype}_dazzler ]] 
        then
            rm -rf ${RAW_FIX_LAFIX_PATH}${ptype}_dazzler
        fi

        mkdir ${FIX_REPMASK_USELAFIX_PATH}_dazzler
        
        ## create a proper dazzler fasta header
        for x in $(seq 1 ${rawblocks})
        do  
            echo "${DACCORD_PATH}/bin/fastaidrename < ${FIX_REPMASK_USELAFIX_PATH}/${RAW_DB%.db}.${x}${RAW_FIX_LAFIX_FILESUFFIX}.fasta | awk '{print \$1}' > ${FIX_REPMASK_USELAFIX_PATH}_dazzler/${RAW_DB%.db}.${x}${RAW_FIX_LAFIX_FILESUFFIX}.fasta"            
        done >> mask_${sID}_createDB_single_${FIX_DB%.db}.${slurmID}.plan        

        # create dazzler db 
        for x in $(seq 1 ${rawblocks})
        do  
        	if [[ $x -eq 1 ]]
        	then 
        		echo "${DAZZLER_PATH}/bin/fasta2DB -v ${FIX_DAZZ_DB%.db}.db ${FIX_REPMASK_USELAFIX_PATH}_dazzler/${RAW_DB%.db}.1${RAW_FIX_LAFIX_FILESUFFIX}.fasta"
        		echo "${DAZZLER_PATH}/bin/DBsplit${SCRUB_DBSPLIT_OPT} ${FIX_DAZZ_DB%.db}.db"                		
        	else
            	echo "${DAZZLER_PATH}/bin/fasta2DB -v ${FIX_DAZZ_DB%.db}.db ${FIX_REPMASK_USELAFIX_PATH}_dazzler/${RAW_DB%.db}.${x}${RAW_FIX_LAFIX_FILESUFFIX}.fasta"
        	fi
        done >> mask_${sID}_createDB_single_${FIX_DB%.db}.${slurmID}.plan
        
        echo "if [[ -d ${FIX_REPMASK_OUTDIR} ]]; then mv ${FIX_REPMASK_OUTDIR} ${FIX_REPMASK_OUTDIR}_$(date '+%Y-%m-%d_%H-%M-%S'); fi && mkdir ${FIX_REPMASK_OUTDIR} && ln -s -r .${FIX_DB%.db}.idx .${FIX_DB%.db}.bps ${FIX_DB%.db}.db .${FIX_DAZZ_DB%.db}.idx .${FIX_DAZZ_DB%.db}.bps ${FIX_DAZZ_DB%.db}.db ${FIX_REPMASK_OUTDIR}" >> mask_${sID}_createDB_single_${FIX_DB%.db}.${slurmID}.plan
        
    	echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mask_${sID}_createDB_single_${FIX_DB%.db}.${slurmID}.version
    	echo "DAZZLER $(git --git-dir=${DAZZLER_SOURCE_PATH}/DAZZ_DB/.git rev-parse --short HEAD)" >> mask_${sID}_createDB_single_${FIX_DB%.db}.${slurmID}.version
    	echo "fastaidrename $(git --git-dir=${DACCORD_SOURCE_PATH}/.git rev-parse --short HEAD)" >> mask_${sID}_createDB_single_${FIX_DB%.db}.${slurmID}.version 
    elif [[ ${currentStep} -eq 2 ]]
    then
        ### clean up plans 
        for x in $(ls mask_${sID}_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        ### find and set DBdust options 
        setDBdustOptions
        ### create DBdust commands 
        for x in $(seq 1 ${fixblocks})
        do 
            echo "cd ${FIX_REPMASK_OUTDIR} && ${MARVEL_PATH}/bin/DBdust${REPMASK_DBDUST_OPT} ${FIX_DB%.db}.${x} && cd ${myCWD}"
            echo "cd ${FIX_REPMASK_OUTDIR} && ${DAZZLER_PATH}/bin/DBdust${REPMASK_DBDUST_OPT} ${FIX_DAZZ_DB%.db}.${x} && cd ${myCWD}"
    	done > mask_${sID}_DBdust_block_${FIX_DB%.db}.${slurmID}.plan
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mask_${sID}_DBdust_block_${FIX_DB%.db}.${slurmID}.version
    	echo "DAZZLER $(git --git-dir=${DAZZLER_SOURCE_PATH}/DAZZ_DB/.git rev-parse --short HEAD)" >> mask_${sID}_DBdust_block_${FIX_DB%.db}.${slurmID}.version        
    elif [[ ${currentStep} -eq 3 ]]
    then 
        ### clean up plans 
        for x in $(ls mask_${sID}_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        ### find and set Catrack options 
        setCatrackOptions
        ### create Catrack command
        echo "cd ${FIX_REPMASK_OUTDIR} && ${MARVEL_PATH}/bin/Catrack${REPMASK_CATRACK_OPT} ${FIX_DB%.db} dust && cp .${FIX_DB%.db}.dust.anno .${FIX_DB%.db}.dust.data ${myCWD}/ && cd ${myCWD}" > mask_${sID}_Catrack_single_${FIX_DB%.db}.${slurmID}.plan
        echo "cd ${FIX_REPMASK_OUTDIR} && ${DAZZLER_PATH}/bin/Catrack${REPMASK_CATRACK_OPT} ${FIX_DAZZ_DB%.db} dust && cp .${FIX_DAZZ_DB%.db}.dust.anno .${FIX_DAZZ_DB%.db}.dust.data ${myCWD}/ && cd ${myCWD}" >> mask_${sID}_Catrack_single_${FIX_DB%.db}.${slurmID}.plan         
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mask_${sID}_Catrack_single_${FIX_DB%.db}.${slurmID}.version
    	echo "DAZZLER $(git --git-dir=${DAZZLER_SOURCE_PATH}/DAZZ_DB/.git rev-parse --short HEAD)" >> mask_${sID}_Catrack_single_${FIX_DB%.db}.${slurmID}.version
    elif [[ ${currentStep} -eq 4 ]]
    then 
        ### clean up plans 
        for x in $(ls mask_${sID}_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done     
        ### find and set datander options 
        setDatanderOptions
        ### create datander commands
        for x in $(seq 1 ${fixblocks})
        do 
            echo "cd ${FIX_REPMASK_OUTDIR} && PATH=${DAZZLER_PATH}/bin:\${PATH} ${DAZZLER_PATH}/bin/datander${REPMASK_DATANDER_OPT} ${FIX_DAZZ_DB%.db}.${x} && cd ${myCWD}"
		done > mask_${sID}_datander_block_${FIX_DB%.db}.${slurmID}.plan
        echo "DAZZLER datander $(git --git-dir=${DAZZLER_SOURCE_PATH}/DAMASKER/.git rev-parse --short HEAD)" > mask_${sID}_datander_block_${FIX_DB%.db}.${slurmID}.version
    elif [[ ${currentStep} -eq 5 ]]
    then 
        ### clean up plans 
        for x in $(ls mask_${sID}_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done     
        ### find and set TANmask options 
        setTANmaskOptions
        ### create TANmask commands
        for x in $(seq 1 ${fixblocks})
        do 
            echo "cd ${FIX_REPMASK_OUTDIR} && ${DAZZLER_PATH}/bin/TANmask${REPMASK_TANMASK_OPT} ${FIX_DAZZ_DB%.db} TAN.${FIX_DAZZ_DB%.db}.${x}.las && cd ${myCWD}" 
    	done > mask_${sID}_TANmask_block_${FIX_DB%.db}.${slurmID}.plan
        echo "DAZZLER TANmask $(git --git-dir=${DAZZLER_SOURCE_PATH}/DAMASKER/.git rev-parse --short HEAD)" > mask_${sID}_TANmask_block_${FIX_DB%.db}.${slurmID}.version
    elif [[ ${currentStep} -eq 6 ]]
    then 
        ### clean up plans 
        for x in $(ls mask_${sID}_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        ### find and set Catrack options
        if [[ -z ${REPMASK_CATRACK_OPT} ]] 
        then
            setCatrackOptions
        fi
        ### create Catrack command
		echo "cd ${FIX_REPMASK_OUTDIR} && ${DAZZLER_PATH}/bin/Catrack${REPMASK_CATRACK_OPT} ${FIX_DAZZ_DB%.db} ${FIX_REPMASK_TANMASK_TRACK} && cp .${FIX_DAZZ_DB%.db}.${FIX_REPMASK_TANMASK_TRACK}.anno .${FIX_DAZZ_DB%.db}.${FIX_REPMASK_TANMASK_TRACK}.data ${myCWD}/ && cd ${myCWD}" > mask_${sID}_Catrack_single_${FIX_DB%.db}.${slurmID}.plan
        echo "cd ${FIX_REPMASK_OUTDIR} && ${LASTOOLS_PATH}/bin/viewmasks ${FIX_DAZZ_DB%.db} ${FIX_REPMASK_TANMASK_TRACK} > ${FIX_DAZZ_DB%.db}.${FIX_REPMASK_TANMASK_TRACK}.txt && cd ${myCWD}" >> mask_${sID}_Catrack_single_${FIX_DB%.db}.${slurmID}.plan
      	echo "cd ${FIX_REPMASK_OUTDIR} && ${MARVEL_PATH}/bin/txt2track -m ${FIX_DB%.db} ${FIX_DAZZ_DB%.db}.${FIX_REPMASK_TANMASK_TRACK}.txt ${FIX_REPMASK_TANMASK_TRACK} && cp .${FIX_DB%.db}.${FIX_REPMASK_TANMASK_TRACK}.a2 .${FIX_DB%.db}.${FIX_REPMASK_TANMASK_TRACK}.d2 ${myCWD}/ && cd ${myCWD}" >> mask_${sID}_Catrack_single_${FIX_DB%.db}.${slurmID}.plan
      	echo "cd ${FIX_REPMASK_OUTDIR} && ${MARVEL_PATH}/bin/TKcombine ${FIX_DB%.db} ${FIX_REPMASK_TANMASK_TRACK}_dust ${FIX_REPMASK_TANMASK_TRACK} dust && cp .${FIX_DB%.db}.${FIX_REPMASK_TANMASK_TRACK}_dust.a2 .${FIX_DB%.db}.${FIX_REPMASK_TANMASK_TRACK}_dust.d2 ${myCWD}/ && cd ${myCWD}" >> mask_${sID}_Catrack_single_${FIX_DB%.db}.${slurmID}.plan 
        
        echo "DAZZLER Catrack $(git --git-dir=${DAZZLER_SOURCE_PATH}/DAZZ_DB/.git rev-parse --short HEAD)" > mask_${sID}_Catrack_single_${FIX_DB%.db}.${slurmID}.version
        echo "LASTOOLS viewmasks $(git --git-dir=${LASTOOLS_SOURCE_PATH}/.git rev-parse --short HEAD)" >> mask_${sID}_Catrack_single_${FIX_DB%.db}.${slurmID}.version    
        echo "DAMAR txt2track $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" >> mask_${sID}_Catrack_single_${FIX_DB%.db}.${slurmID}.version
        echo "DAMAR TKcombine $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" >> mask_${sID}_Catrack_single_${FIX_DB%.db}.${slurmID}.version        
    elif [[ ${currentStep} -eq 7 ]]
    then
        for x in $(ls mask_${sID}_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        ### find and set daligner options 
        setDalignerOptions

		## create job directories before daligner runs
		for x in $(seq 1 ${fixblocks})
		do
			if [[ -d ${FIX_REPMASK_OUTDIR}/mask_${x}_B${FIX_REPMASK_BLOCKCMP[0]}C${FIX_REPMASK_LAREPEAT_COV[0]} ]]
			then
				mv ${FIX_REPMASK_OUTDIR}/mask_${x}_B${FIX_REPMASK_BLOCKCMP[0]}C${FIX_REPMASK_LAREPEAT_COV[0]} ${RAW_REPAMSK_OUTDIR}/mask_${x}_B${FIX_REPMASK_BLOCKCMP[0]}C${FIX_REPMASK_LAREPEAT_COV[0]}_$(date '+%Y-%m-%d_%H-%M-%S')	
			fi
			mkdir -p ${FIX_REPMASK_OUTDIR}/mask_${x}_B${FIX_REPMASK_BLOCKCMP[0]}C${FIX_REPMASK_LAREPEAT_COV[0]}	
		done		

        bcmp=${FIX_REPMASK_BLOCKCMP[0]}

        ### create daligner commands
        n=${bcmp}
        for x in $(seq 1 ${fixblocks})
        do
            if [[ $(echo "$x%${bcmp}" | bc) -eq 1 || ${bcmp} -eq 1 ]]
            then 
              n=${bcmp}
            fi 
            if [[ -n ${FIX_REPMASK_REPEATTRACK} ]]
            then
                REP="-m${FIX_REPMASK_REPEATTRACK}"
            fi
            if [[ -n ${FIX_REPMASK_DALIGNER_NUMACTL} && ${FIX_REPMASK_DALIGNER_NUMACTL} -gt 0 ]] && [[ "x${SLURM_NUMACTL}" == "x" || ${SLURM_NUMACTL} -eq 0 ]]
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
            echo -n "cd ${FIX_REPMASK_OUTDIR} && PATH=${DAZZLER_PATH}/bin:\${PATH} ${NUMACTL}${DAZZLER_PATH}/bin/daligner${REPMASK_DALIGNER_OPT} ${REP} ${FIX_DAZZ_DB%.db}.${x}"
            for y in $(seq ${x} $((${x}+${n}-1)))
            do
                if [[ ${y} -gt ${fixblocks} ]]
                then
                    break
                fi
                echo -n " ${FIX_DAZZ_DB%.db}.${y}"
            done 
            
			for y in $(seq ${x} $((${x}+${n}-1)))
            do
                if [[ ${y} -gt ${fixblocks} ]]
                then
                    break
                fi
                echo -n " && mv ${FIX_DAZZ_DB%.db}.${x}.${FIX_DAZZ_DB%.db}.${y}.las mask_${x}_B${FIX_REPMASK_BLOCKCMP[0]}C${FIX_REPMASK_LAREPEAT_COV[0]}"
            done 
            
            n=$((${n}-1))

            echo " && cd ${myCWD}"
    	done > mask_${sID}_daligner_block_${FIX_DB%.db}.${slurmID}.plan
        echo "DAZZLER daligner $(git --git-dir=${DAZZLER_SOURCE_PATH}/DALIGNER/.git rev-parse --short HEAD)" > mask_${sID}_daligner_block_${FIX_DB%.db}.${slurmID}.version
    elif [[ ${currentStep} -eq 8 ]]
    then
        ### clean up plans 
        for x in $(ls mask_${sID}_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 

        ### create LAmerge commands 
        for x in $(seq 1 ${fixblocks})
        do 
        	echo "cd ${FIX_REPMASK_OUTDIR} && ${MARVEL_PATH}/bin/LAmerge -n 32 ${FIX_DB%.db} ${FIX_DAZZ_DB%.db}.${x}.maskB${FIX_REPMASK_BLOCKCMP[0]}C${FIX_REPMASK_LAREPEAT_COV[0]}.las mask_${x}_B${FIX_REPMASK_BLOCKCMP[0]}C${FIX_REPMASK_LAREPEAT_COV[0]} && cd ${myCWD}"
    	done > mask_${sID}_LAmerge_block_${FIX_DB%.db}.${slurmID}.plan      
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mask_${sID}_LAmerge_block_${FIX_DB%.db}.${slurmID}.version  
    elif [[ ${currentStep} -eq 9 ]]
    then 
        ### clean up plans 
        for x in $(ls mask_${sID}_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        ### find and set LArepeat options 
        setLArepeatOptions 0
        ### create LArepeat commands
        for x in $(seq 1 ${fixblocks})
        do 
        	echo "cd ${FIX_REPMASK_OUTDIR} && ${MARVEL_PATH}/bin/LArepeat${REPMASK_LAREPEAT_OPT} -b ${x} ${FIX_DB%.db} ${FIX_DAZZ_DB%.db}.${x}.maskB${FIX_REPMASK_BLOCKCMP[0]}C${FIX_REPMASK_LAREPEAT_COV[0]}.las && cd ${myCWD}/" 
            echo "cd ${FIX_REPMASK_OUTDIR} && ln -s -f ${FIX_DAZZ_DB%.db}.${x}.maskB${FIX_REPMASK_BLOCKCMP[0]}C${FIX_REPMASK_LAREPEAT_COV[0]}.las ${FIX_DAZZ_DB%.db}.${x}.maskB${FIX_REPMASK_BLOCKCMP[0]}C${FIX_REPMASK_LAREPEAT_COV[0]}.${x}.las && ${DAZZLER_PATH}/bin/REPmask -v -c${FIX_REPMASK_LAREPEAT_COV[0]} -n${FIX_REPMASK_REPEATTRACK} ${FIX_DAZZ_DB%.db} ${FIX_DAZZ_DB%.db}.${x}.maskB${FIX_REPMASK_BLOCKCMP[0]}C${FIX_REPMASK_LAREPEAT_COV[0]}.${x}.las && unlink ${FIX_DAZZ_DB%.db}.${x}.maskB${FIX_REPMASK_BLOCKCMP[0]}C${FIX_REPMASK_LAREPEAT_COV[0]}.${x}.las && cd ${myCWD}/"
		done > mask_${sID}_LArepeat_block_${FIX_DB%.db}.${slurmID}.plan
		echo "MARVEL LArepeat $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mask_${sID}_LArepeat_block_${FIX_DB%.db}.${slurmID}.version
        echo "DAZZLER REPmask $(git --git-dir=${DAZZLER_SOURCE_PATH}/DAMASKER/.git rev-parse --short HEAD)" >> mask_${sID}_LArepeat_block_${FIX_DB%.db}.${slurmID}.version
    elif [[ ${currentStep} -eq 10 ]]
    then 
        ### clean up plans 
        for x in $(ls mask_${sID}_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        ### find and set TKmerge options 
        setTKmergeOptions
        setLArepeatOptions 0
        ### create TKmerge commands
        echo "cd ${FIX_REPMASK_OUTDIR} && ${MARVEL_PATH}/bin/TKmerge${REPMASK_TKMERGE_OPT} ${FIX_DB%.db} ${FIX_REPMASK_REPEATTRACK} && cp .${FIX_DB%.db}.${FIX_REPMASK_REPEATTRACK}.a2 .${FIX_DB%.db}.${FIX_REPMASK_REPEATTRACK}.d2 ${myCWD}/ && cd ${myCWD}/" > mask_${sID}_TKmerge_single_${FIX_DB%.db}.${slurmID}.plan
        echo "cd ${FIX_REPMASK_OUTDIR} && ${DAZZLER_PATH}/bin/Catrack${REPMASK_TKMERGE_OPT} -f -v ${FIX_DAZZ_DB%.db} ${FIX_REPMASK_REPEATTRACK} && cp .${FIX_DAZZ_DB%.db}.${FIX_REPMASK_REPEATTRACK}.anno .${FIX_DAZZ_DB%.db}.${FIX_REPMASK_REPEATTRACK}.data ${myCWD}/ && cd ${myCWD}/" >> mask_${sID}_TKmerge_single_${FIX_DB%.db}.${slurmID}.plan
        echo "MARVEL TKmerge $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mask_${sID}_TKmerge_single_${FIX_DB%.db}.${slurmID}.version
        echo "DAZZLER Catrack $(git --git-dir=${DAZZLER_SOURCE_PATH}/DAZZ_DB/.git rev-parse --short HEAD)" >> mask_${sID}_TKmerge_single_${FIX_DB%.db}.${slurmID}.version    
    elif [[ ${currentStep} -eq 11 && ${#FIX_REPMASK_BLOCKCMP[*]} -eq 2 && ${#FIX_REPMASK_LAREPEAT_COV[*]} -eq 2 ]]
    then
        ### clean up plans 
        for x in $(ls mask_${sID}_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        ### find and set daligner options 
        setDalignerOptions

        setLArepeatOptions 0
        bcmp=${FIX_REPMASK_BLOCKCMP[1]}

		## create job directories before daligner runs
		for x in $(seq 1 ${fixblocks})
		do
			if [[ -d ${FIX_REPMASK_OUTDIR}/mask_${x}_B${FIX_REPMASK_BLOCKCMP[1]}C${FIX_REPMASK_LAREPEAT_COV[1]} ]]
			then
				mv ${FIX_REPMASK_OUTDIR}/mask_${x}_B${FIX_REPMASK_BLOCKCMP[1]}C${FIX_REPMASK_LAREPEAT_COV[1]} ${FIX_REPMASK_OUTDIR}/mask_${x}_B${FIX_REPMASK_BLOCKCMP[1]}C${FIX_REPMASK_LAREPEAT_COV[1]}_$(date '+%Y-%m-%d_%H-%M-%S')	
			fi
			mkdir -p ${RAW_REPAMSK_OUTDIR}/mask_${x}_B${FIX_REPMASK_BLOCKCMP[1]}C${FIX_REPMASK_LAREPEAT_COV[1]}	
		done

        ### create daligner commands
        n=${bcmp}
        for x in $(seq 1 ${fixblocks})
        do
            if [[ $(echo "$x%${bcmp}" | bc) -eq 1 || ${bcmp} -eq 1 ]]
            then 
              n=$((${bcmp}))
            fi 
            if [[ -n ${FIX_REPMASK_REPEATTRACK} ]]
            then
                REP="-m${FIX_REPMASK_REPEATTRACK}"
            fi
            if [[ -n ${FIX_REPMASK_DALIGNER_NUMACTL} && ${FIX_REPMASK_DALIGNER_NUMACTL} -gt 0 ]] && [[ "x${SLURM_NUMACTL}" == "x" || ${SLURM_NUMACTL} -eq 0 ]]
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

           if [[ "x${DALIGNER_VERSION}" == "x2" ]]
			then
				echo -n "cd ${FIX_REPMASK_OUTDIR} && PATH=${DAZZLER_PATH}/bin:\${PATH} ${NUMACTL}${DAZZLER_PATH}/bin/daligner${REPMASK_DALIGNER_OPT} ${REP} ${FIX_DAZZ_DB%.db}.${x} ${FIX_DAZZ_DB%.db}.@${x}"
			else
				echo -n "cd ${FIX_REPMASK_OUTDIR} && PATH=${DAZZLER_PATH}/bin:\${PATH} ${NUMACTL}${DAZZLER_PATH}/bin/daligner${REPMASK_DALIGNER_OPT} ${REP} ${FIX_DAZZ_DB%.db}.${x}"
			fi			
			
            for y in $(seq ${x} $((${x}+${n}-1)))
            do
                if [[ ${y} -gt ${fixblocks} ]]
                then
                	y=$((y-1))
                    break
                fi
                if [[ "x${DALIGNER_VERSION}" != "x2" ]]
				then
					echo -n " ${FIX_DAZZ_DB%.db}.${y}"
				fi			
                                
            done 
            
            if [[ "x${DALIGNER_VERSION}" == "x2" ]]
			then
				echo -n "-${y} && mv"
			else
				echo -n " && mv"
			fi			
            
            for y in $(seq ${x} $((${x}+${n}-1)))
            do
                if [[ ${y} -gt ${fixblocks} ]]
                then
                    break
                fi
                echo -n " ${FIX_DAZZ_DB%.db}.${x}.${FIX_DAZZ_DB%.db}.${y}.las"
            done
            echo -n " mask_${x}_B${FIX_REPMASK_BLOCKCMP[1]}C${FIX_REPMASK_LAREPEAT_COV[1]}"
            
            
			if [[ -z "${FIX_REPMASK_DALIGNER_ASYMMETRIC}" || ${FIX_REPMASK_DALIGNER_ASYMMETRIC} -ne 0 ]]
			then
				
				for y in $(seq $((x+1)) $((x+n-1)))
            	do
                	if [[ ${y} -gt ${fixblocks} ]]
                	then
                    	break
                	fi
                	echo -n " && mv ${FIX_DAZZ_DB%.db}.${y}.${FIX_DAZZ_DB%.db}.${x}.las mask_${y}_B${FIX_REPMASK_BLOCKCMP[1]}C${FIX_REPMASK_LAREPEAT_COV[1]}"
            	done
        	fi
 
            echo " && cd ${myCWD}"
            n=$((${n}-1))
    	done > mask_${sID}_daligner_block_${FIX_DB%.db}.${slurmID}.plan 
        echo "DAZZLER daligner $(git --git-dir=${DAZZLER_SOURCE_PATH}/DALIGNER/.git rev-parse --short HEAD)" > mask_${sID}_daligner_block_${FIX_DB%.db}.${slurmID}.version
    elif [[ ${currentStep} -eq 12 && ${#FIX_REPMASK_BLOCKCMP[*]} -eq 2 && ${#FIX_REPMASK_LAREPEAT_COV[*]} -eq 2 ]]
    then
        ### clean up plans 
        for x in $(ls mask_${sID}_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 

        ### create LAmerge commands 
        for x in $(seq 1 ${fixblocks})
        do 
        	 echo "cd ${FIX_REPMASK_OUTDIR} && ${MARVEL_PATH}/bin/LAmerge -n 32 ${FIX_DB%.db} ${FIX_DAZZ_DB%.db}.${x}.maskB${FIX_REPMASK_BLOCKCMP[1]}C${FIX_REPMASK_LAREPEAT_COV[1]}.las mask_${x}_B${FIX_REPMASK_BLOCKCMP[1]}C${FIX_REPMASK_LAREPEAT_COV[1]} && cd ${myCWD}"
    	done > mask_${sID}_LAmerge_block_${FIX_DB%.db}.${slurmID}.plan
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mask_${sID}_LAmerge_block_${FIX_DB%.db}.${slurmID}.version
    elif [[ ${currentStep} -eq 13 && ${#FIX_REPMASK_BLOCKCMP[*]} -eq 2 && ${#FIX_REPMASK_LAREPEAT_COV[*]} -eq 2 ]]
    then 
        ### clean up plans 
        for x in $(ls mask_${sID}_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        ### find and set LArepeat options 
        setLArepeatOptions 1
        ### create LArepeat commands
        for x in $(seq 1 ${fixblocks})
        do 
            echo "cd ${FIX_REPMASK_OUTDIR} && ${MARVEL_PATH}/bin/LArepeat${REPMASK_LAREPEAT_OPT} -b ${x} ${FIX_DB%.db} ${FIX_DAZZ_DB%.db}.${x}.maskB${FIX_REPMASK_BLOCKCMP[1]}C${FIX_REPMASK_LAREPEAT_COV[1]}.las && cd ${myCWD}/" 
            echo "cd ${FIX_REPMASK_OUTDIR} && ln -s -f ${FIX_DAZZ_DB%.db}.${x}.maskB${FIX_REPMASK_BLOCKCMP[1]}C${FIX_REPMASK_LAREPEAT_COV[1]}.las ${FIX_DAZZ_DB%.db}.${x}.maskB${FIX_REPMASK_BLOCKCMP[1]}C${FIX_REPMASK_LAREPEAT_COV[1]}.${x}.las && ${DAZZLER_PATH}/bin/REPmask -v -c${FIX_REPMASK_LAREPEAT_COV[1]} -n${FIX_REPMASK_REPEATTRACK} ${FIX_DAZZ_DB%.db} ${FIX_DAZZ_DB%.db}.${x}.maskB${FIX_REPMASK_BLOCKCMP[1]}C${FIX_REPMASK_LAREPEAT_COV[1]}.${x}.las && unlink ${FIX_DAZZ_DB%.db}.${x}.maskB${FIX_REPMASK_BLOCKCMP[1]}C${FIX_REPMASK_LAREPEAT_COV[1]}.${x}.las && cd ${myCWD}/"
    	done > mask_${sID}_LArepeat_block_${FIX_DB%.db}.${slurmID}.plan
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mask_${sID}_LArepeat_block_${FIX_DB%.db}.${slurmID}.version
    elif [[ ${currentStep} -eq 14 && ${#FIX_REPMASK_BLOCKCMP[*]} -eq 2 && ${#FIX_REPMASK_LAREPEAT_COV[*]} -eq 2 ]]
    then 
        ### clean up plans 
        for x in $(ls mask_${sID}_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        ### find and set TKmerge options 
        setTKmergeOptions
        setLArepeatOptions 1
        ### create TKmerge commands
        echo "cd ${FIX_REPMASK_OUTDIR} && ${MARVEL_PATH}/bin/TKmerge${REPMASK_TKMERGE_OPT} ${FIX_DB%.db} ${FIX_REPMASK_REPEATTRACK} && cp .${FIX_DB%.db}.${FIX_REPMASK_REPEATTRACK}.a2 .${FIX_DB%.db}.${FIX_REPMASK_REPEATTRACK}.d2 ${myCWD}/ && cd ${myCWD}/" > mask_${sID}_TKmerge_single_${FIX_DB%.db}.${slurmID}.plan
        echo "cd ${FIX_REPMASK_OUTDIR} && ${DAZZLER_PATH}/bin/Catrack${REPMASK_TKMERGE_OPT} -f -v ${FIX_DAZZ_DB%.db} ${FIX_REPMASK_REPEATTRACK} && cp .${FIX_DAZZ_DB%.db}.${FIX_REPMASK_REPEATTRACK}.anno .${FIX_DAZZ_DB%.db}.${FIX_REPMASK_REPEATTRACK}.data ${myCWD}/ && cd ${myCWD}/" >> mask_${sID}_TKmerge_single_${FIX_DB%.db}.${slurmID}.plan
        echo "MARVEL TKmerge $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mask_${sID}_TKmerge_single_${FIX_DB%.db}.${slurmID}.version
        echo "DAZZLER Catrack $(git --git-dir=${DAZZLER_SOURCE_PATH}/DAZZ_DB/.git rev-parse --short HEAD)" >> mask_${sID}_TKmerge_single_${FIX_DB%.db}.${slurmID}.version
    else 
        (>&2 echo "step ${currentStep} in FIX_REPMASK_TYPE ${FIX_REPMASK_TYPE} not supported")
        (>&2 echo "valid steps are: ${myTypes[${FIX_REPMASK_TYPE}]}")
        exit 1        
    fi    
else
    echo "unknown FIX_REPMASK_TYPE ${FIX_REPMASK_TYPE}"
    (>&2 echo "supported types")
    x=0; while [ $x -lt ${#myTypes[*]} ]; do (>&2 echo "${myTypes[${x}]}"); done
    exit 1
fi

exit 0