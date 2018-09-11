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
        REPMASK_TANMASK_OPT="${REPMASK_TANMASK_OPT} -l ${FIX_REPMASK_TANMASK_MINLEN}"
    fi
    if [[ -n ${FIX_REPMASK_TANMASK_TRACK} ]]
    then
        REPMASK_TANMASK_OPT="${REPMASK_TANMASK_OPT} -m ${FIX_REPMASK_TANMASK_TRACK}"
    fi
}

function setDaligerOptions()
{
    REPMASK_DALIGNER_OPT=""
    if [[ -n ${FIX_REPMASK_DALIGNER_IDENTITY_OVLS} && ${FIX_REPMASK_DALIGNER_IDENTITY_OVLS} -gt 0 ]]
    then
        REPMASK_DALIGNER_OPT="${REPMASK_DALIGNER_OPT} -I"
    fi
    if [[ -n ${FIX_REPMASK_DALIGNER_KMER} && ${FIX_REPMASK_DALIGNER_KMER} -gt 0 ]]
    then
        REPMASK_DALIGNER_OPT="${REPMASK_DALIGNER_OPT} -k ${FIX_REPMASK_DALIGNER_KMER}"
    fi
    if [[ -n ${FIX_REPMASK_DALIGNER_ERR} ]]
    then
        REPMASK_DALIGNER_OPT="${REPMASK_DALIGNER_OPT} -e ${FIX_REPMASK_DALIGNER_ERR}"
    fi
    if [[ -n ${FIX_REPMASK_DALIGNER_BIAS} && ${FIX_REPMASK_DALIGNER_BIAS} -eq 1 ]]
    then
        REPMASK_DALIGNER_OPT="${REPMASK_DALIGNER_OPT} -b"
    fi
    if [[ -n ${FIX_REPMASK_DALIGNER_RUNID} ]]
    then
        REPMASK_DALIGNER_OPT="${REPMASK_DALIGNER_OPT} -r ${FIX_REPMASK_DALIGNER_RUNID}"
    fi
    if [[ -n ${FIX_REPMASK_DALIGNER_OLEN} ]]
    then
        REPMASK_DALIGNER_OPT="${REPMASK_DALIGNER_OPT} -l ${FIX_REPMASK_DALIGNER_OLEN}"
    fi    
    if [[ -n ${FIX_REPMASK_DALIGNER_MEM} && ${FIX_REPMASK_DALIGNER_MEM} -gt 0 ]]
    then
        REPMASK_DALIGNER_OPT="${REPMASK_DALIGNER_OPT} -M ${FIX_REPMASK_DALIGNER_MEM}"
    fi    
    if [[ -n ${FIX_REPMASK_DALIGNER_HITS} ]]
    then
        REPMASK_DALIGNER_OPT="${REPMASK_DALIGNER_OPT} -h ${FIX_REPMASK_DALIGNER_HITS}"
    fi        
    if [[ -n ${FIX_REPMASK_DALIGNER_T} ]]
    then
        REPMASK_DALIGNER_OPT="${REPMASK_DALIGNER_OPT} -t ${FIX_REPMASK_DALIGNER_T}"
    fi  
    if [[ -n ${FIX_REPMASK_DALIGNER_MASK} ]]
    then
        for x in ${FIX_REPMASK_DALIGNER_MASK}
        do 
            REPMASK_DALIGNER_OPT="${REPMASK_DALIGNER_OPT} -m ${x}"
        done
    fi
    if [[ -n ${FIX_REPMASK_DALIGNER_TRACESPACE} && ${FIX_REPMASK_DALIGNER_TRACESPACE} -gt 0 ]]
    then
        REPMASK_DALIGNER_OPT="${REPMASK_DALIGNER_OPT} -s ${FIX_REPMASK_DALIGNER_TRACESPACE}"
    fi
    if [[ -n ${THREADS_daligner} ]]
    then 
        REPMASK_DALIGNER_OPT="${REPMASK_DALIGNER_OPT} -j ${THREADS_daligner}"
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
        REPMASK_DATANDER_OPT="${REPMASK_DATANDER_OPT} -j ${FIX_REPMASK_DATANDER_THREADS}"
    fi
    if [[ -n ${FIX_REPMASK_DATANDER_MINLEN} ]]
    then
        REPMASK_DATANDER_OPT="${REPMASK_DATANDER_OPT} -l ${FIX_REPMASK_DATANDER_MINLEN}"
    fi
    if [[ -n ${FIX_REPMASK_DATANDER_FOLDER} ]]
    then
        REPMASK_DATANDER_OPT="${REPMASK_DATANDER_OPT} -o ${FIX_REPMASK_DATANDER_FOLDER}"
    else
        FIX_REPMASK_DATANDER_FOLDER="tan"
        REPMASK_DATANDER_OPT="${REPMASK_DATANDER_OPT} -o ${FIX_REPMASK_DATANDER_FOLDER}"
    fi
}

function setDBsplitOptions()
{
    REPMASK_DBSPLIT_OPT=""

    FIX_REPMASK_DBSPLIT_S=$(grep size ${RAW_DB%.db}.db | awk '{print $3}')

    if [[ -z ${FIX_REPMASK_DBSPLIT_S} ]]
    then
        (>&2 echo "database $db has not been partitioned. Run DBsplit first!")
        exit 1
    fi

    REPMASK_DBSPLIT_OPT="${REPMASK_DBSPLIT_OPT} -s${FIX_REPMASK_DBSPLIT_S}"
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

myTypes=("1-createFIX_DB, 2-DBdust, 3-Catrack, 4-datander, 5-TANmask, 6-Catrack, 7-daligner, 8-LAmerge, 9-LArepeat, 10-TKmerge, 11-daligner, 12-LAmerge, 13-LArepeat, 14-TKmerge")
# type_0 - steps: 1-createFIX_DB, 2-DBdust, 3-Catrack, 4-datander, 5-TANmask, 6-Catrack, 7-daligner, 8-LAmerge, 9-LArepeat, 10-TKmerge, 11-daligner, 12-LAmerge, 13-LArepeat, 14-TKmerge
if [[ ${FIX_REPMASK_TYPE} -eq 0 ]]
then
	if [[ ${currentStep} -eq 1 ]]
    then
    	### clean up plans 
        for x in $(ls mask_01_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 

        if [[ -f ${FIX_DB%.db}.db ]]; 
        then 
            echo "p3_s1 rm existing DB ${FIX_DB%.db}.db"
            ${MARVEL_PATH}/bin/DBrm ${FIX_DB%.db}.db            
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
        done > mask_01_createDB_single_${FIX_DB%.db}.${slurmID}.plan        

        ### convert corrected reads into a proper dazzler read format and create a valid dazzler db, that where read IDS do exactly map to the corressponding marvel db

        if [[ -f ${FIX_DAZZ_DB%.db}.db ]]; 
        then 
            echo "p3_s1 rm existing DAZZ DB ${FIX_DAZZ_DB%.db}.db"
            ${MARVEL_PATH}/bin/DBrm ${FIX_DAZZ_DB%.db}.db            
        fi

        if [[ -d ${RAW_FIX_LAFIX_PATH}${ptype}_dazzler ]] 
        then
            rm -rf ${RAW_FIX_LAFIX_PATH}${ptype}_dazzler
        fi

        mkdir ${FIX_REPMASK_USELAFIX_PATH}_dazzler
        
        ## create a proper dazzler fasta header
        for x in $(seq 1 ${rawblocks})
        do  
            echo "${DACCORD_PATH}/bin/fastaidrename < ${FIX_REPMASK_USELAFIX_PATH}/${RAW_DB%.db}.${x}${RAW_FIX_LAFIX_FILESUFFIX}.fasta | awk '{print \$1}' > ${FIX_REPMASK_USELAFIX_PATH}_dazzler/${RAW_DB%.db}.${x}${RAW_FIX_LAFIX_FILESUFFIX}.fasta"            
        done >> mask_01_createDB_single_${FIX_DB%.db}.${slurmID}.plan        

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
        done >> mask_01_createDB_single_${FIX_DB%.db}.${slurmID}.plan
    	echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mask_01_createDB_single_${FIX_DB%.db}.${slurmID}.version
    	echo "DAZZLER $(git --git-dir=${DAZZLER_SOURCE_PATH}/DAZZ_DB/.git rev-parse --short HEAD)" >> mask_01_createDB_single_${FIX_DB%.db}.${slurmID}.version
    	echo "fastaidrename $(git --git-dir=${DACCORD_SOURCE_PATH}/.git rev-parse --short HEAD)" >> mask_01_createDB_single_${FIX_DB%.db}.${slurmID}.version 
    elif [[ ${currentStep} -eq 2 ]]
    then
        ### clean up plans 
        for x in $(ls mask_02_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        ### find and set DBdust options 
        setDBdustOptions
        ### create DBdust commands 
        for x in $(seq 1 ${fixblocks})
        do 
            echo "${MARVEL_PATH}/bin/DBdust${REPMASK_DBDUST_OPT} ${FIX_DB%.db}.${x}"
            echo "${DAZZLER_PATH}/bin/DBdust${REPMASK_DBDUST_OPT} ${FIX_DAZZ_DB%.db}.${x}"
    	done > mask_02_DBdust_block_${FIX_DB%.db}.${slurmID}.plan
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mask_02_DBdust_block_${FIX_DB%.db}.${slurmID}.version
    	echo "DAZZLER $(git --git-dir=${DAZZLER_SOURCE_PATH}/DAZZ_DB/.git rev-parse --short HEAD)" >> mask_02_DBdust_block_${FIX_DB%.db}.${slurmID}.version        
    elif [[ ${currentStep} -eq 3 ]]
    then 
        ### clean up plans 
        for x in $(ls mask_03_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        ### find and set Catrack options 
        setCatrackOptions
        ### create Catrack command
        echo "${MARVEL_PATH}/bin/Catrack${REPMASK_CATRACK_OPT} ${FIX_DB%.db} dust" > mask_03_Catrack_single_${FIX_DB%.db}.${slurmID}.plan
        echo "${DAZZLER_PATH}/bin/Catrack${REPMASK_CATRACK_OPT} ${FIX_DAZZ_DB%.db} dust" >> mask_03_Catrack_single_${FIX_DB%.db}.${slurmID}.plan         
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mask_03_Catrack_single_${FIX_DB%.db}.${slurmID}.version
    	echo "DAZZLER $(git --git-dir=${DAZZLER_SOURCE_PATH}/DAZZ_DB/.git rev-parse --short HEAD)" >> mask_02_DBdust_block_${FIX_DB%.db}.${slurmID}.version
    elif [[ ${currentStep} -eq 4 ]]
    then 
        ### clean up plans 
        for x in $(ls mask_04_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done     
        ### find and set datander options 
        setDatanderOptions
        ### create datander commands
        for x in $(seq 1 ${fixblocks})
        do 
            echo "${MARVEL_PATH}/bin/datander${REPMASK_DATANDER_OPT} ${FIX_DB%.db}.${x}"
		done > mask_04_datander_block_${FIX_DB%.db}.${slurmID}.plan
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mask_04_datander_block_${FIX_DB%.db}.${slurmID}.version
    elif [[ ${currentStep} -eq 5 ]]
    then 
        ### clean up plans 
        for x in $(ls mask_05_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done     
        ### find and set TANmask options 
        setTANmaskOptions
        ### create TANmask commands
        for x in $(seq 1 ${fixblocks})
        do 
            echo "${MARVEL_PATH}/bin/TANmask${REPMASK_TANMASK_OPT} ${FIX_DB%.db} ${FIX_REPMASK_DATANDER_FOLDER}/${FIX_DB%.db}.${x}.${FIX_DB%.db}.${x}.las" 
    	done > mask_05_TANmask_block_${FIX_DB%.db}.${slurmID}.plan
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mask_05_TANmask_block_${FIX_DB%.db}.${slurmID}.version
    elif [[ ${currentStep} -eq 6 ]]
    then 
        ### clean up plans 
        for x in $(ls mask_06_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        ### find and set Catrack options
        if [[ -z ${REPMASK_CATRACK_OPT} ]] 
        then
            setCatrackOptions
        fi
        ### create Catrack command
        echo "${MARVEL_PATH}/bin/Catrack${REPMASK_CATRACK_OPT} ${FIX_DB%.db} ${FIX_REPMASK_TANMASK_TRACK}" > mask_06_Catrack_single_${FIX_DB%.db}.${slurmID}.plan
        echo "${MARVEL_PATH}/bin/TKcombine ${FIX_DB%.db} ${FIX_REPMASK_TANMASK_TRACK}_dust ${FIX_REPMASK_TANMASK_TRACK} dust" >> mask_06_Catrack_single_${FIX_DB%.db}.${slurmID}.plan
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mask_06_Catrack_single_${FIX_DB%.db}.${slurmID}.version    
    elif [[ ${currentStep} -eq 7 ]]
    then
        for x in $(ls mask_07_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        ### find and set daligner options 
        setDaligerOptions

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
            if [[ -n ${FIX_REPMASK_DALIGNER_NUMACTL} && ${FIX_REPMASK_DALIGNER_NUMACTL} -gt 0 ]]
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
            echo -n "${NUMACTL}${MARVEL_PATH}/bin/daligner${REPMASK_DALIGNER_OPT} ${REP} ${FIX_DB%.db}.${x}"
            for y in $(seq ${x} $((${x}+${n}-1)))
            do
                if [[ ${y} -gt ${fixblocks} ]]
                then
                    break
                fi
                echo -n " ${FIX_DB%.db}.${y}"
            done 
            n=$((${n}-1))

            echo ""
    	done > mask_07_daligner_block_${FIX_DB%.db}.${slurmID}.plan
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mask_07_daligner_block_${FIX_DB%.db}.${slurmID}.version
    elif [[ ${currentStep} -eq 8 ]]
    then
        ### clean up plans 
        for x in $(ls mask_08_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 

        ### create LAmerge commands 
        for x in $(seq 1 ${fixblocks})
        do 
            echo "${MARVEL_PATH}/bin/LAmerge -n 32 ${FIX_DB%.db} ${FIX_DB%.db}.${x}.mask1.las $(getSubDirName ${FIX_REPMASK_DALIGNER_RUNID} ${x})"
    	done > mask_08_LAmerge_block_${FIX_DB%.db}.${slurmID}.plan      
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mask_08_LAmerge_block_${FIX_DB%.db}.${slurmID}.version  
    elif [[ ${currentStep} -eq 9 ]]
    then 
        ### clean up plans 
        for x in $(ls mask_09_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        ### find and set LArepeat options 
        setLArepeatOptions 0
        ### create LArepeat commands
        for x in $(seq 1 ${fixblocks})
        do 
            echo "${MARVEL_PATH}/bin/LArepeat${REPMASK_LAREPEAT_OPT} -b ${x} ${FIX_DB%.db} ${FIX_DB%.db}.${x}.mask1.las"
		done > mask_09_LArepeat_block_${FIX_DB%.db}.${slurmID}.plan
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mask_09_LArepeat_block_${FIX_DB%.db}.${slurmID}.version
    elif [[ ${currentStep} -eq 10 ]]
    then 
        ### clean up plans 
        for x in $(ls mask_10_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        ### find and set TKmerge options 
        setTKmergeOptions
        setLArepeatOptions 0
        ### create TKmerge commands
        echo "${MARVEL_PATH}/bin/TKmerge${REPMASK_TKMERGE_OPT} ${FIX_DB%.db} ${FIX_REPMASK_REPEATTRACK}" > mask_10_TKmerge_single_${FIX_DB%.db}.${slurmID}.plan
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mask_10_TKmerge_single_${FIX_DB%.db}.${slurmID}.version    
    elif [[ ${currentStep} -eq 11 && ${#FIX_REPMASK_BLOCKCMP[*]} -eq 2 && ${#FIX_REPMASK_LAREPEAT_COV[*]} -eq 2 ]]
    then
        ### clean up plans 
        for x in $(ls mask_11_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        ### find and set daligner options 
        setDaligerOptions

        setLArepeatOptions 0
        bcmp=${FIX_REPMASK_BLOCKCMP[1]}

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
            if [[ -n ${FIX_REPMASK_DALIGNER_NUMACTL} && ${FIX_REPMASK_DALIGNER_NUMACTL} -gt 0 ]]
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

            echo -n "${NUMACTL}${MARVEL_PATH}/bin/daligner${REPMASK_DALIGNER_OPT} ${REP} ${FIX_DB%.db}.${x}"
            for y in $(seq ${x} $((${x}+${n}-1)))
            do
                if [[ ${y} -gt ${fixblocks} ]]
                then
                    break
                fi
                echo -n " ${FIX_DB%.db}.${y}"
            done 
            n=$((${n}-1))

            echo ""
    	done > mask_11_daligner_block_${FIX_DB %.db}.${slurmID}.plan 
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mask_11_daligner_block_${FIX_DB%.db}.${slurmID}.version
    elif [[ ${currentStep} -eq 12 && ${#FIX_REPMASK_BLOCKCMP[*]} -eq 2 && ${#FIX_REPMASK_LAREPEAT_COV[*]} -eq 2 ]]
    then
        ### clean up plans 
        for x in $(ls mask_12_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 

        ### create LAmerge commands 
        for x in $(seq 1 ${fixblocks})
        do 
            echo "${MARVEL_PATH}/bin/LAmerge -n 32 ${FIX_DB%.db} ${FIX_DB%.db}.${x}.mask2.las $(getSubDirName ${FIX_REPMASK_DALIGNER_RUNID} ${x})"
    	done > mask_12_LAmerge_block_${FIX_DB %.db}.${slurmID}.plan
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mask_12_LAmerge_block_${FIX_DB%.db}.${slurmID}.version
    elif [[ ${currentStep} -eq 13 && ${#FIX_REPMASK_BLOCKCMP[*]} -eq 2 && ${#FIX_REPMASK_LAREPEAT_COV[*]} -eq 2 ]]
    then 
        ### clean up plans 
        for x in $(ls mask_13_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        ### find and set LArepeat options 
        setLArepeatOptions 1
        ### create LArepeat commands
        for x in $(seq 1 ${fixblocks})
        do 
            echo "${MARVEL_PATH}/bin/LArepeat${REPMASK_LAREPEAT_OPT} -b ${x} ${FIX_DB%.db} ${FIX_DB%.db}.${x}.mask2.las"
    	done > mask_13_LArepeat_block_${FIX_DB%.db}.${slurmID}.plan
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mask_13_LArepeat_block_${FIX_DB%.db}.${slurmID}.version
    elif [[ ${currentStep} -eq 14 && ${#FIX_REPMASK_BLOCKCMP[*]} -eq 2 && ${#FIX_REPMASK_LAREPEAT_COV[*]} -eq 2 ]]
    then 
        ### clean up plans 
        for x in $(ls mask_14_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        ### find and set TKmerge options 
        setTKmergeOptions
        setLArepeatOptions 1
        ### create TKmerge commands
        echo "${MARVEL_PATH}/bin/TKmerge${REPMASK_TKMERGE_OPT} ${FIX_DB%.db} ${FIX_REPMASK_REPEATTRACK}" > mask_14_TKmerge_single_${FIX_DB%.db}.${slurmID}.plan
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mask_14_TKmerge_single_${FIX_DB%.db}.${slurmID}.version    
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