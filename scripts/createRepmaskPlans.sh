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

if [[ ! -n "${RAW_REPMASK_TYPE}" ]]
then 
    (>&2 echo "cannot create repmask jobs if varibale RAW_REPMASK_TYPE is not set.")
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

nblocks=$(getNumOfDbBlocks)

function setDBdustOptions()
{
    REPMASK_DBDUST_OPT=""
    if [[ -n ${RAW_REPMASK_DBDUST_BIAS} && ${RAW_REPMASK_DBDUST_BIAS} -ge 1 ]]
    then
        REPMASK_DBDUST_OPT="${REPMASK_DBDUST_OPT} -b"
    fi
}

function setCatrackOptions()
{
    REPMASK_CATRACK_OPT=""
    if [[ -n ${RAW_REPMASK_CATRACK_VERBOSE} && ${RAW_REPMASK_CATRACK_VERBOSE} -ge 1 ]]
    then
        REPMASK_CATRACK_OPT="${REPMASK_CATRACK_OPT} -v"
    fi
    if [[ -n ${RAW_REPMASK_CATRACK_DELETE} && ${RAW_REPMASK_CATRACK_DELETE} -ge 1 ]]
    then
        REPMASK_CATRACK_OPT="${REPMASK_CATRACK_OPT} -d"
    fi
    if [[ -n ${RAW_REPMASK_CATRACK_OVERWRITE} && ${RAW_REPMASK_CATRACK_OVERWRITE} -ge 1 ]]
    then
        REPMASK_CATRACK_OPT="${REPMASK_CATRACK_OPT} -f"
    fi
}

function setTANmaskOptions()
{
    REPMASK_TANMASK_OPT=""
    if [[ -n ${RAW_REPMASK_TANMASK_VERBOSE} && ${RAW_REPMASK_TANMASK_VERBOSE} -ge 1 ]]
    then
        REPMASK_TANMASK_OPT="${REPMASK_TANMASK_OPT} -v"
    fi
    if [[ -n ${RAW_REPMASK_TANMASK_MINLEN} && ${RAW_REPMASK_TANMASK_MINLEN} -ge 1 ]]
    then
        REPMASK_TANMASK_OPT="${REPMASK_TANMASK_OPT} -l${RAW_REPMASK_TANMASK_MINLEN}"
    fi
    if [[ -n ${RAW_REPMASK_TANMASK_TRACK} ]]
    then
        REPMASK_TANMASK_OPT="${REPMASK_TANMASK_OPT} -n${RAW_REPMASK_TANMASK_TRACK}"
    fi
}

function setDaligerOptions()
{
    REPMASK_DALIGNER_OPT=""
    if [[ -n ${RAW_REPMASK_DALIGNER_IDENTITY_OVLS} && ${RAW_REPMASK_DALIGNER_IDENTITY_OVLS} -gt 0 ]]
    then
        REPMASK_DALIGNER_OPT="${REPMASK_DALIGNER_OPT} -I"
    fi
    if [[ -n ${RAW_REPMASK_DALIGNER_KMER} && ${RAW_REPMASK_DALIGNER_KMER} -gt 0 ]]
    then
        REPMASK_DALIGNER_OPT="${REPMASK_DALIGNER_OPT} -k${RAW_REPMASK_DALIGNER_KMER}"
    fi
    if [[ -n ${RAW_REPMASK_DALIGNER_ERR} ]]
    then
        REPMASK_DALIGNER_OPT="${REPMASK_DALIGNER_OPT} -e${RAW_REPMASK_DALIGNER_ERR}"
    fi
    if [[ -n ${RAW_REPMASK_DALIGNER_BIAS} && ${RAW_REPMASK_DALIGNER_BIAS} -eq 1 ]]
    then
        REPMASK_DALIGNER_OPT="${REPMASK_DALIGNER_OPT} -b"
    fi
    if [[ -n ${RAW_REPMASK_DALIGNER_BRIDGE} && ${RAW_REPMASK_DALIGNER_BRIDGE} -eq 1 ]]
    then
        REPMASK_DALIGNER_OPT="${REPMASK_DALIGNER_OPT} -B"
    fi
    if [[ -n ${RAW_REPMASK_DALIGNER_OLEN} ]]
    then
        REPMASK_DALIGNER_OPT="${REPMASK_DALIGNER_OPT} -l${RAW_REPMASK_DALIGNER_OLEN}"
    fi    
    if [[ -n ${RAW_REPMASK_DALIGNER_MEM} && ${RAW_REPMASK_DALIGNER_MEM} -gt 0 ]]
    then
        REPMASK_DALIGNER_OPT="${REPMASK_DALIGNER_OPT} -M${RAW_REPMASK_DALIGNER_MEM}"
    fi    
    if [[ -n ${RAW_REPMASK_DALIGNER_HITS} ]]
    then
        REPMASK_DALIGNER_OPT="${REPMASK_DALIGNER_OPT} -h${RAW_REPMASK_DALIGNER_HITS}"
    fi        
    if [[ -n ${RAW_REPMASK_DALIGNER_T} ]]
    then
        REPMASK_DALIGNER_OPT="${REPMASK_DALIGNER_OPT} -t${RAW_REPMASK_DALIGNER_T}"
    fi  
    if [[ -n ${RAW_REPMASK_DALIGNER_MASK} ]]
    then
        for x in ${RAW_REPMASK_DALIGNER_MASK}
        do 
            REPMASK_DALIGNER_OPT="${REPMASK_DALIGNER_OPT} -m${x}"
        done
    fi
    if [[ -n ${RAW_REPMASK_DALIGNER_TRACESPACE} && ${RAW_REPMASK_DALIGNER_TRACESPACE} -gt 0 ]]
    then
        REPMASK_DALIGNER_OPT="${REPMASK_DALIGNER_OPT} -s${RAW_REPMASK_DALIGNER_TRACESPACE}"
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
    if [[ -n ${RAW_REPMASK_LAREPEAT_LOW} ]]
    then
        REPMASK_LAREPEAT_OPT="${REPMASK_LAREPEAT_OPT} -l ${RAW_REPMASK_LAREPEAT_LOW}"
    fi
    if [[ -n ${RAW_REPMASK_LAREPEAT_HGH} ]]
    then
        REPMASK_LAREPEAT_OPT="${REPMASK_LAREPEAT_OPT} -h ${RAW_REPMASK_LAREPEAT_HGH}"
    fi
    if [[ -n ${RAW_REPMASK_LAREPEAT_OLEN} ]]
    then
        REPMASK_LAREPEAT_OPT="${REPMASK_LAREPEAT_OPT} -o ${RAW_REPMASK_LAREPEAT_OLEN}"
    fi
    if [[ -n ${RAW_REPMASK_LAREPEAT_REPEATTRACK} ]]
    then
        if [[ -z ${idx} ]]
        then
          RAW_REPMASK_REPEATTRACK=${RAW_REPMASK_LAREPEAT_REPEATTRACK}
          REPMASK_LAREPEAT_OPT="${REPMASK_LAREPEAT_OPT} -t ${RAW_REPMASK_REPEATTRACK}"
        else 
            if [[ ${#RAW_REPMASK_LAREPEAT_COV[*]} -lt ${idx} ]]
            then 
                (>&2 echo "RAW_REPMASK_LAREPEAT_COV has lower the ${idx} elements")
                exit 1
            elif [[ ${#RAW_REPMASK_BLOCKCMP[*]} -lt ${idx} ]]
            then 
                (>&2 echo "RAW_REPMASK_BLOCKCMP has lower the ${idx} elements")
                exit 1
            fi
            RAW_REPMASK_REPEATTRACK=${RAW_REPMASK_LAREPEAT_REPEATTRACK}_B${RAW_REPMASK_BLOCKCMP[${idx}]}C${RAW_REPMASK_LAREPEAT_COV[${idx}]}
            REPMASK_LAREPEAT_OPT="${REPMASK_LAREPEAT_OPT} -t ${RAW_REPMASK_REPEATTRACK}"
        fi
    fi
    if [[ -n ${RAW_REPMASK_LAREPEAT_COV} ]]
    then
        if [[ -z ${idx} ]]
        then
            REPMASK_LAREPEAT_OPT="${REPMASK_LAREPEAT_OPT} -c ${RAW_REPMASK_LAREPEAT_COV}"
        else
            if [[ ${#RAW_REPMASK_LAREPEAT_COV[*]} -lt ${idx} ]]
            then 
                (>&2 echo "RAW_REPMASK_LAREPEAT_COV has lower the ${idx} elements")
                exit 1
            fi
            REPMASK_LAREPEAT_OPT="${REPMASK_LAREPEAT_OPT} -c ${RAW_REPMASK_LAREPEAT_COV[${idx}]}"
        fi 
    fi
}

function setTKmergeOptions()
{
    REPMASK_TKMERGE_OPT=""
    if [[ -n ${RAW_REPMASK_TKMERGE_DELETE} && ${RAW_REPMASK_TKMERGE_DELETE} -ge 1 ]]
    then
        REPMASK_TKMERGE_OPT="${REPMASK_TKMERGE_OPT} -d"
    fi
    if [ ! -n ${RAW_REPMASK_LAREPEAT_REPEATTRACK} ] ### fall back to default value!!!
    then
        RAW_REPMASK_LAREPEAT_REPEATTRACK="repeats"
    fi
}

function setDatanderOptions()
{
    ### find and set datander options 
    REPMASK_DATANDER_OPT=""
    if [[ -n ${RAW_REPMASK_DATANDER_THREADS} ]]
    then
        REPMASK_DATANDER_OPT="${REPMASK_DATANDER_OPT} -T${RAW_REPMASK_DATANDER_THREADS}"
    fi
    if [[ -n ${RAW_REPMASK_DATANDER_MINLEN} ]]
    then
        REPMASK_DATANDER_OPT="${REPMASK_DATANDER_OPT} -l${RAW_REPMASK_DATANDER_MINLEN}"
    fi
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

if [[ -z "${RAW_REPMASK_OUTDIR}" ]]
then
	RAW_REPMASK_OUTDIR=repmask	
fi

myCWD=$(pwd)

if [[ ${currentStep} -lt 10 ]]
then 
	sID=0${currentStep}
else
	sID=${currentStep}
fi

myTypes=("01_createSubdir, 02_DBdust, 03_Catrack, 04_datander, 05_TANmask, 06_Catrack, 07_daligner, 08_LAmerge, 09_LArepeat, 10_TKmerge, 11-daligner, 12-LAmerge, 13-LArepeat, 14-TKmerge")
# type_0 - stepsp[1-14}: 01_createSubdir, 02_DBdust, 03_Catrack, 04_datander, 05_TANmask, 06_Catrack, 07_daligner, 08_LAmerge, 09_LArepeat, 10_TKmerge, 11-daligner, 12-LAmerge, 13-LArepeat, 14-TKmerge
if [[ ${RAW_REPMASK_TYPE} -eq 0 ]]
then
	if [[ ${currentStep} -eq 1 ]]
    then
        ### clean up plans 
        for x in $(ls mask_${sID}_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        echo "if [[ -d ${RAW_REPMASK_OUTDIR} ]]; then mv ${RAW_REPMASK_OUTDIR} ${RAW_REPMASK_OUTDIR}_$(date '+%Y-%m-%d_%H-%M-%S'); fi && mkdir ${RAW_REPMASK_OUTDIR} && ln -s -r .${RAW_DB%.db}.idx .${RAW_DB%.db}.bps ${RAW_DB%.db}.db .${RAW_DAZZ_DB%.db}.idx .${RAW_DAZZ_DB%.db}.bps ${RAW_DAZZ_DB%.db}.db ${RAW_REPMASK_OUTDIR}" > mask_${sID}_createSubdir_single_${RAW_DB%.db}.${slurmID}.plan
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mask_${sID}_createSubdir_single_${RAW_DB%.db}.${slurmID}.version         
    elif [[ ${currentStep} -eq 2 ]]
    then
        ### clean up plans 
        for x in $(ls mask_${sID}_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        ### find and set DBdust options 
        setDBdustOptions
        
        ### create DBdust commands 
        for x in $(seq 1 ${nblocks})
        do 
            echo "cd ${RAW_REPMASK_OUTDIR} && ${MARVEL_PATH}/bin/DBdust${REPMASK_DBDUST_OPT} ${RAW_DB%.db}.${x} && cd ${myCWD}"
            echo "cd ${RAW_REPMASK_OUTDIR} && ${DAZZLER_PATH}/bin/DBdust${REPMASK_DBDUST_OPT} ${RAW_DAZZ_DB%.db}.${x} && cd ${myCWD}"
    	done > mask_${sID}_DBdust_block_${RAW_DB%.db}.${slurmID}.plan
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mask_${sID}_DBdust_block_${RAW_DB%.db}.${slurmID}.version
        echo "DAZZLER $(git --git-dir=${DAZZLER_SOURCE_PATH}/DAZZ_DB/.git rev-parse --short HEAD)" >> mask_${sID}_DBdust_block_${RAW_DB%.db}.${slurmID}.version
    elif [[ ${currentStep} -eq 3 ]]
    then 
        ### clean up plans 
        for x in $(ls mask_03_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        
        ### find and set Catrack options 
        setCatrackOptions
        ### create Catrack command
        echo "cd ${RAW_REPMASK_OUTDIR} && ${MARVEL_PATH}/bin/Catrack${REPMASK_CATRACK_OPT} ${RAW_DB%.db} dust && cp .${RAW_DB%.db}.dust.anno .${RAW_DB%.db}.dust.data ${myCWD}/ && cd ${myCWD}" > mask_03_Catrack_single_${RAW_DB%.db}.${slurmID}.plan
        echo "cd ${RAW_REPMASK_OUTDIR} && ${DAZZLER_PATH}/bin/Catrack${REPMASK_CATRACK_OPT} ${RAW_DAZZ_DB%.db} dust && cp .${RAW_DAZZ_DB%.db}.dust.anno .${RAW_DAZZ_DB%.db}.dust.data ${myCWD}/ && cd ${myCWD}" >> mask_03_Catrack_single_${RAW_DB%.db}.${slurmID}.plan
                 
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mask_03_Catrack_single_${RAW_DB%.db}.${slurmID}.version
        echo "DAZZLER $(git --git-dir=${DAZZLER_SOURCE_PATH}/DAZZ_DB/.git rev-parse --short HEAD)" >> mask_03_Catrack_single_${RAW_DB%.db}.${slurmID}.version
    elif [[ ${currentStep} -eq 4 ]]
    then 
        ### clean up plans 
        for x in $(ls mask_04_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done     
        ### find and set datander options 
        setDatanderOptions
        
        ### create datander commands
        for x in $(seq 1 ${nblocks})
        do 
            echo "cd ${RAW_REPMASK_OUTDIR} && PATH=${DAZZLER_PATH}/bin:\${PATH} ${DAZZLER_PATH}/bin/datander${REPMASK_DATANDER_OPT} ${RAW_DAZZ_DB%.db}.${x} && cd ${myCWD}"
    	done > mask_04_datander_block_${RAW_DB%.db}.${slurmID}.plan
        echo "DAZZLER datander $(git --git-dir=${DAZZLER_SOURCE_PATH}/DAMASKER/.git rev-parse --short HEAD)" > mask_04_datander_block_${RAW_DB%.db}.${slurmID}.version
    elif [[ ${currentStep} -eq 5 ]]
    then 
        ### clean up plans 
        for x in $(ls mask_05_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
             
        ### find and set TANmask options         
        setTANmaskOptions
        ### create TANmask commands
        for x in $(seq 1 ${nblocks})
        do 
            echo "cd ${RAW_REPMASK_OUTDIR} && ${DAZZLER_PATH}/bin/TANmask${REPMASK_TANMASK_OPT} ${RAW_DAZZ_DB%.db} TAN.${RAW_DAZZ_DB%.db}.${x}.las && cd ${myCWD}" 
    	done > mask_05_TANmask_block_${RAW_DB%.db}.${slurmID}.plan
        echo "DAZZLER TANmask $(git --git-dir=${DAZZLER_SOURCE_PATH}/DAMASKER/.git rev-parse --short HEAD)" > mask_05_TANmask_block_${RAW_DB%.db}.${slurmID}.version
    elif [[ ${currentStep} -eq 6 ]]
    then 
        ### clean up plans 
        for x in $(ls mask_06_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        
        ### find and set Catrack options
        if [[ -z ${REPMASK_CATRACK_OPT} ]] 
        then
            setCatrackOptions
        fi
        ### create Catrack command
        echo "cd ${RAW_REPMASK_OUTDIR} && ${DAZZLER_PATH}/bin/Catrack${REPMASK_CATRACK_OPT} ${RAW_DAZZ_DB%.db} ${RAW_REPMASK_TANMASK_TRACK} && cp .${RAW_DAZZ_DB%.db}.${RAW_REPMASK_TANMASK_TRACK}.anno .${RAW_DAZZ_DB%.db}.${RAW_REPMASK_TANMASK_TRACK}.data ${myCWD}/ && cd ${myCWD}" > mask_06_Catrack_single_${RAW_DB%.db}.${slurmID}.plan
        echo "cd ${RAW_REPMASK_OUTDIR} && ${DAZZLER_PATH}/bin/DBdump -r -m${RAW_REPMASK_TANMASK_TRACK} ${RAW_DAZZ_DB%.db} | awk '{if (\$1 == \"R\") {read=\$2}; if (\$1 == \"T0\" && \$2 > 0) {for (i = 3; i < 3+2*\$2; i+=2) print read-1\" \"\$i\" \"\$(i+1)} }' > ${RAW_DAZZ_DB%.db}.${RAW_REPMASK_TANMASK_TRACK}.txt && cd ${myCWD}" >> mask_06_Catrack_single_${RAW_DB%.db}.${slurmID}.plan
      	echo "cd ${RAW_REPMASK_OUTDIR} && ${MARVEL_PATH}/bin/txt2track -m ${RAW_DB%.db} ${RAW_DAZZ_DB%.db}.${RAW_REPMASK_TANMASK_TRACK}.txt ${RAW_REPMASK_TANMASK_TRACK} && cp .${RAW_DB%.db}.${RAW_REPMASK_TANMASK_TRACK}.a2 .${RAW_DB%.db}.${RAW_REPMASK_TANMASK_TRACK}.d2 ${myCWD}/ && cd ${myCWD}" >> mask_06_Catrack_single_${RAW_DB%.db}.${slurmID}.plan
      	echo "cd ${RAW_REPMASK_OUTDIR} && ${MARVEL_PATH}/bin/TKcombine ${RAW_DB%.db} ${RAW_REPMASK_TANMASK_TRACK}_dust ${RAW_REPMASK_TANMASK_TRACK} dust && cp .${RAW_DB%.db}.${RAW_REPMASK_TANMASK_TRACK}_dust.a2 .${RAW_DB%.db}.${RAW_REPMASK_TANMASK_TRACK}_dust.d2 ${myCWD}/ && cd ${myCWD}" >> mask_06_Catrack_single_${RAW_DB%.db}.${slurmID}.plan 
        
        echo "DAZZLER Catrack $(git --git-dir=${DAZZLER_SOURCE_PATH}/DAZZ_DB/.git rev-parse --short HEAD)" > mask_06_Catrack_single_${RAW_DB%.db}.${slurmID}.version
        echo "LASTOOLS viewmasks $(git --git-dir=${LASTOOLS_SOURCE_PATH}/.git rev-parse --short HEAD)" >> mask_06_Catrack_single_${RAW_DB%.db}.${slurmID}.version    
        echo "DAMAR txt2track $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" >> mask_06_Catrack_single_${RAW_DB%.db}.${slurmID}.version
        echo "DAMAR TKcombine $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" >> mask_06_Catrack_single_${RAW_DB%.db}.${slurmID}.version
    elif [[ ${currentStep} -eq 7 ]]
    then
        for x in $(ls mask_07_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        ### find and set daligner options 
        setDaligerOptions
		
		## create job directories before daligner runs
		for x in $(seq 1 ${nblocks})
		do
			if [[ -d ${RAW_REPMASK_OUTDIR}/mask_${x}_B${RAW_REPMASK_BLOCKCMP[0]}C${RAW_REPMASK_LAREPEAT_COV[0]} ]]
			then
				mv ${RAW_REPMASK_OUTDIR}/mask_${x}_B${RAW_REPMASK_BLOCKCMP[0]}C${RAW_REPMASK_LAREPEAT_COV[0]} ${RAW_REPMASK_OUTDIR}/mask_${x}_B${RAW_REPMASK_BLOCKCMP[0]}C${RAW_REPMASK_LAREPEAT_COV[0]}_$(date '+%Y-%m-%d_%H-%M-%S')	
			fi
			mkdir -p ${RAW_REPMASK_OUTDIR}/mask_${x}_B${RAW_REPMASK_BLOCKCMP[0]}C${RAW_REPMASK_LAREPEAT_COV[0]}	
		done

        bcmp=${RAW_REPMASK_BLOCKCMP[0]}
		
        ### create daligner commands
        n=${bcmp}
        for x in $(seq 1 ${nblocks})
        do
            if [[ $(echo "$x%${bcmp}" | bc) -eq 1 || ${bcmp} -eq 1 ]]
            then 
              n=${bcmp}
            fi 
            if [[ -n ${RAW_REPMASK_REPEATTRACK} ]]
            then
                REP="-m${RAW_REPMASK_REPEATTRACK}"
            fi
            if [[ -n ${RAW_REPMASK_DALIGNER_NUMACTL} && ${RAW_REPMASK_DALIGNER_NUMACTL} -gt 0 ]] && [[ "x${SLURM_NUMACTL}" == "x" || ${SLURM_NUMACTL} -eq 0 ]]
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
            echo -n "cd ${RAW_REPMASK_OUTDIR} && PATH=${DAZZLER_PATH}/bin:\${PATH} ${NUMACTL}${DAZZLER_PATH}/bin/daligner${REPMASK_DALIGNER_OPT} ${REP} ${RAW_DAZZ_DB%.db}.${x}"
            for y in $(seq ${x} $((${x}+${n}-1)))
            do
                if [[ ${y} -gt ${nblocks} ]]
                then
                    break
                fi
                echo -n " ${RAW_DAZZ_DB%.db}.${y}"
            done 

			for y in $(seq ${x} $((${x}+${n}-1)))
            do
                if [[ ${y} -gt ${nblocks} ]]
                then
                    break
                fi
                echo -n " && mv ${RAW_DAZZ_DB%.db}.${x}.${RAW_DAZZ_DB%.db}.${y}.las mask_${x}_B${RAW_REPMASK_BLOCKCMP[0]}C${RAW_REPMASK_LAREPEAT_COV[0]}"
            done 
            
            n=$((${n}-1))

            echo " && cd ${myCWD}"
   		done > mask_07_daligner_block_${RAW_DB%.db}.${slurmID}.plan
        echo "DAZZLER daligner $(git --git-dir=${DAZZLER_SOURCE_PATH}/DALIGNER/.git rev-parse --short HEAD)" > mask_07_daligner_block_${RAW_DB%.db}.${slurmID}.version
    elif [[ ${currentStep} -eq 8 ]]
    then
        ### clean up plans 
        for x in $(ls mask_08_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
		
        ### create LAmerge commands 
        for x in $(seq 1 ${nblocks})
        do 
            echo "cd ${RAW_REPMASK_OUTDIR} && ${MARVEL_PATH}/bin/LAmerge -n 32 ${RAW_DB%.db} ${RAW_DAZZ_DB%.db}.${x}.maskB${RAW_REPMASK_BLOCKCMP[0]}C${RAW_REPMASK_LAREPEAT_COV[0]}.las mask_${x}_B${RAW_REPMASK_BLOCKCMP[0]}C${RAW_REPMASK_LAREPEAT_COV[0]} && cd ${myCWD}"            
    	done > mask_08_LAmerge_block_${RAW_DB%.db}.${slurmID}.plan
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mask_08_LAmerge_block_${RAW_DB%.db}.${slurmID}.version  
    elif [[ ${currentStep} -eq 9 ]]
    then 
        ### clean up plans 
        for x in $(ls mask_09_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        
        ### find and set LArepeat options 
        setLArepeatOptions 0
        ### create LArepeat commands
        for x in $(seq 1 ${nblocks})
        do 
            echo "cd ${RAW_REPMASK_OUTDIR} && ${MARVEL_PATH}/bin/LArepeat${REPMASK_LAREPEAT_OPT} -b ${x} ${RAW_DB%.db} ${RAW_DAZZ_DB%.db}.${x}.maskB${RAW_REPMASK_BLOCKCMP[0]}C${RAW_REPMASK_LAREPEAT_COV[0]}.las && cd ${myCWD}/" 
            echo "cd ${RAW_REPMASK_OUTDIR} && ln -s -f ${RAW_DAZZ_DB%.db}.${x}.maskB${RAW_REPMASK_BLOCKCMP[0]}C${RAW_REPMASK_LAREPEAT_COV[0]}.las ${RAW_DAZZ_DB%.db}.${x}.maskB${RAW_REPMASK_BLOCKCMP[0]}C${RAW_REPMASK_LAREPEAT_COV[0]}.${x}.las && ${DAZZLER_PATH}/bin/REPmask -v -c${RAW_REPMASK_LAREPEAT_COV[0]} -n${RAW_REPMASK_REPEATTRACK} ${RAW_DAZZ_DB%.db} ${RAW_DAZZ_DB%.db}.${x}.maskB${RAW_REPMASK_BLOCKCMP[0]}C${RAW_REPMASK_LAREPEAT_COV[0]}.${x}.las && unlink ${RAW_DAZZ_DB%.db}.${x}.maskB${RAW_REPMASK_BLOCKCMP[0]}C${RAW_REPMASK_LAREPEAT_COV[0]}.${x}.las && cd ${myCWD}/"
    	done > mask_09_LArepeat_block_${RAW_DB%.db}.${slurmID}.plan
        echo "MARVEL LArepeat $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mask_09_LArepeat_block_${RAW_DB%.db}.${slurmID}.version
        echo "DAZZLER REPmask $(git --git-dir=${DAZZLER_SOURCE_PATH}/DAMASKER/.git rev-parse --short HEAD)" >> mask_09_LArepeat_block_${RAW_DB%.db}.${slurmID}.version
    elif [[ ${currentStep} -eq 10 ]]
    then 
        ### clean up plans 
        for x in $(ls mask_10_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        
        ### find and set TKmerge options 
        setTKmergeOptions
        setLArepeatOptions 0
        ### create TKmerge commands
        echo "cd ${RAW_REPMASK_OUTDIR} && ${MARVEL_PATH}/bin/TKmerge${REPMASK_TKMERGE_OPT} ${RAW_DB%.db} ${RAW_REPMASK_REPEATTRACK} && cp .${RAW_DB%.db}.${RAW_REPMASK_REPEATTRACK}.a2 .${RAW_DB%.db}.${RAW_REPMASK_REPEATTRACK}.d2 ${myCWD}/ && cd ${myCWD}/" > mask_10_TKmerge_single_${RAW_DB%.db}.${slurmID}.plan
        echo "cd ${RAW_REPMASK_OUTDIR} && ${DAZZLER_PATH}/bin/Catrack${REPMASK_TKMERGE_OPT} -f -v ${RAW_DAZZ_DB%.db} ${RAW_REPMASK_REPEATTRACK} && cp .${RAW_DAZZ_DB%.db}.${RAW_REPMASK_REPEATTRACK}.anno .${RAW_DAZZ_DB%.db}.${RAW_REPMASK_REPEATTRACK}.data ${myCWD}/ && cd ${myCWD}/" >> mask_10_TKmerge_single_${RAW_DB%.db}.${slurmID}.plan
        echo "MARVEL TKmerge $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mask_10_TKmerge_single_${RAW_DB%.db}.${slurmID}.version
        echo "DAZZLER Catrack $(git --git-dir=${DAZZLER_SOURCE_PATH}/DAZZ_DB/.git rev-parse --short HEAD)" >> mask_10_TKmerge_single_${RAW_DB%.db}.${slurmID}.version    
    elif [[ ${currentStep} -eq 11 && ${#RAW_REPMASK_BLOCKCMP[*]} -eq 2 && ${#RAW_REPMASK_LAREPEAT_COV[*]} -eq 2 ]]
    then
        ### clean up plans 
        for x in $(ls mask_11_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        ### find and set daligner options 
        setDaligerOptions
		
        setLArepeatOptions 0
        bcmp=${RAW_REPMASK_BLOCKCMP[1]}
			
		## create job directories before daligner runs
		for x in $(seq 1 ${nblocks})
		do
			if [[ -d ${RAW_REPMASK_OUTDIR}/mask_${x}_B${RAW_REPMASK_BLOCKCMP[1]}C${RAW_REPMASK_LAREPEAT_COV[1]} ]]
			then
				mv ${RAW_REPMASK_OUTDIR}/mask_${x}_B${RAW_REPMASK_BLOCKCMP[1]}C${RAW_REPMASK_LAREPEAT_COV[1]} ${RAW_REPMASK_OUTDIR}/mask_${x}_B${RAW_REPMASK_BLOCKCMP[1]}C${RAW_REPMASK_LAREPEAT_COV[1]}_$(date '+%Y-%m-%d_%H-%M-%S')	
			fi
			mkdir -p ${RAW_REPMASK_OUTDIR}/mask_${x}_B${RAW_REPMASK_BLOCKCMP[1]}C${RAW_REPMASK_LAREPEAT_COV[1]}	
		done
		

        ### create daligner commands
        n=${bcmp}
        for x in $(seq 1 ${nblocks})
        do
            if [[ $(echo "$x%${bcmp}" | bc) -eq 1 || ${bcmp} -eq 1 ]]
            then 
              n=$((${bcmp}))
            fi 
            if [[ -n ${RAW_REPMASK_REPEATTRACK} ]]
            then
                REP="-m${RAW_REPMASK_REPEATTRACK}"
            fi
            if [[ -n ${RAW_REPMASK_DALIGNER_NUMACTL} && ${RAW_REPMASK_DALIGNER_NUMACTL} -gt 0 ]] && [[ "x${SLURM_NUMACTL}" == "x" || ${SLURM_NUMACTL} -eq 0 ]]
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
				echo -n "cd ${RAW_REPMASK_OUTDIR} && PATH=${DAZZLER_PATH}/bin:\${PATH} ${NUMACTL}${DAZZLER_PATH}/bin/daligner${REPMASK_DALIGNER_OPT} ${REP} ${RAW_DAZZ_DB%.db}.${x} ${RAW_DAZZ_DB%.db}.@${x}"
			else
				echo -n "cd ${RAW_REPMASK_OUTDIR} && PATH=${DAZZLER_PATH}/bin:\${PATH} ${NUMACTL}${DAZZLER_PATH}/bin/daligner${REPMASK_DALIGNER_OPT} ${REP} ${RAW_DAZZ_DB%.db}.${x}"
			fi			
			
            for y in $(seq ${x} $((${x}+${n}-1)))
            do
                if [[ ${y} -gt ${nblocks} ]]
                then
                	y=$((y-1))
                    break
                fi
                if [[ "x${DALIGNER_VERSION}" != "x2" ]]
				then
					echo -n " ${RAW_DAZZ_DB%.db}.${y}"
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
                if [[ ${y} -gt ${nblocks} ]]
                then
                    break
                fi
                echo -n " ${RAW_DAZZ_DB%.db}.${x}.${RAW_DAZZ_DB%.db}.${y}.las"
            done
            echo -n " mask_${x}_B${RAW_REPMASK_BLOCKCMP[1]}C${RAW_REPMASK_LAREPEAT_COV[1]}"
            
            
			if [[ -z "${RAW_REPMASK_DALIGNER_ASYMMETRIC}" || ${RAW_REPMASK_DALIGNER_ASYMMETRIC} -ne 0 ]]
			then
				
				for y in $(seq $((x+1)) $((x+n-1)))
            	do
                	if [[ ${y} -gt ${nblocks} ]]
                	then
                    	break
                	fi
                	echo -n " && mv ${RAW_DAZZ_DB%.db}.${y}.${RAW_DAZZ_DB%.db}.${x}.las mask_${y}_B${RAW_REPMASK_BLOCKCMP[1]}C${RAW_REPMASK_LAREPEAT_COV[1]}"
            	done
        	fi
 
            echo " && cd ${myCWD}"
            n=$((${n}-1))
    	done > mask_11_daligner_block_${RAW_DB%.db}.${slurmID}.plan 
        echo "DAZZLER daligner $(git --git-dir=${DAZZLER_SOURCE_PATH}/DALIGNER/.git rev-parse --short HEAD)" > mask_11_daligner_block_${RAW_DB%.db}.${slurmID}.version
    elif [[ ${currentStep} -eq 12 && ${#RAW_REPMASK_BLOCKCMP[*]} -eq 2 && ${#RAW_REPMASK_LAREPEAT_COV[*]} -eq 2 ]]
    then
        ### clean up plans 
        for x in $(ls mask_12_*_*_${RAW_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 

        ### create LAmerge commands 
        for x in $(seq 1 ${nblocks})
        do 
            echo "cd ${RAW_REPMASK_OUTDIR} && ${MARVEL_PATH}/bin/LAmerge -n 32 ${RAW_DB%.db} ${RAW_DAZZ_DB%.db}.${x}.maskB${RAW_REPMASK_BLOCKCMP[1]}C${RAW_REPMASK_LAREPEAT_COV[1]}.las mask_${x}_B${RAW_REPMASK_BLOCKCMP[1]}C${RAW_REPMASK_LAREPEAT_COV[1]} && cd ${myCWD}"            
    	done > mask_12_LAmerge_block_${RAW_DB%.db}.${slurmID}.plan
        echo "MARVEL LAmerge  $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mask_11_LAmerge_block_${RAW_DB%.db}.${slurmID}.version
    elif [[ ${currentStep} -eq 13 && ${#RAW_REPMASK_BLOCKCMP[*]} -eq 2 && ${#RAW_REPMASK_LAREPEAT_COV[*]} -eq 2 ]]
    then 
        ### clean up plans 
        for x in $(ls mask_13_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        
        ### find and set LArepeat options 
        setLArepeatOptions 1
        ### create LArepeat commands
        for x in $(seq 1 ${nblocks})
        do 
            echo "cd ${RAW_REPMASK_OUTDIR} && ${MARVEL_PATH}/bin/LArepeat${REPMASK_LAREPEAT_OPT} -b ${x} ${RAW_DB%.db} ${RAW_DAZZ_DB%.db}.${x}.maskB${RAW_REPMASK_BLOCKCMP[1]}C${RAW_REPMASK_LAREPEAT_COV[1]}.las && cd ${myCWD}/" 
            echo "cd ${RAW_REPMASK_OUTDIR} && ln -s -f ${RAW_DAZZ_DB%.db}.${x}.maskB${RAW_REPMASK_BLOCKCMP[1]}C${RAW_REPMASK_LAREPEAT_COV[1]}.las ${RAW_DAZZ_DB%.db}.${x}.maskB${RAW_REPMASK_BLOCKCMP[1]}C${RAW_REPMASK_LAREPEAT_COV[1]}.${x}.las && ${DAZZLER_PATH}/bin/REPmask -v -c${RAW_REPMASK_LAREPEAT_COV[1]} -n${RAW_REPMASK_REPEATTRACK} ${RAW_DAZZ_DB%.db} ${RAW_DAZZ_DB%.db}.${x}.maskB${RAW_REPMASK_BLOCKCMP[1]}C${RAW_REPMASK_LAREPEAT_COV[1]}.${x}.las && unlink ${RAW_DAZZ_DB%.db}.${x}.maskB${RAW_REPMASK_BLOCKCMP[1]}C${RAW_REPMASK_LAREPEAT_COV[1]}.${x}.las && cd ${myCWD}/"
    	done > mask_13_LArepeat_block_${RAW_DB%.db}.${slurmID}.plan
        echo "MARVEL LArepeat $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mask_13_LArepeat_block_${RAW_DB%.db}.${slurmID}.version
        echo "DAZZLER REPmask $(git --git-dir=${DAZZLER_SOURCE_PATH}/DAMASKER/.git rev-parse --short HEAD)" >> mask_13_LArepeat_block_${RAW_DB%.db}.${slurmID}.version
    elif [[ ${currentStep} -eq 14 && ${#RAW_REPMASK_BLOCKCMP[*]} -eq 2 && ${#RAW_REPMASK_LAREPEAT_COV[*]} -eq 2 ]]
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
        echo "cd ${RAW_REPMASK_OUTDIR} && ${MARVEL_PATH}/bin/TKmerge${REPMASK_TKMERGE_OPT} ${RAW_DB%.db} ${RAW_REPMASK_REPEATTRACK} && cp .${RAW_DB%.db}.${RAW_REPMASK_REPEATTRACK}.a2 .${RAW_DB%.db}.${RAW_REPMASK_REPEATTRACK}.d2 ${myCWD}/ && cd ${myCWD}/" > mask_14_TKmerge_single_${RAW_DB%.db}.${slurmID}.plan
        echo "cd ${RAW_REPMASK_OUTDIR} && ${DAZZLER_PATH}/bin/Catrack${REPMASK_TKMERGE_OPT} -f -v ${RAW_DAZZ_DB%.db} ${RAW_REPMASK_REPEATTRACK} && cp .${RAW_DAZZ_DB%.db}.${RAW_REPMASK_REPEATTRACK}.anno .${RAW_DAZZ_DB%.db}.${RAW_REPMASK_REPEATTRACK}.data ${myCWD}/ && cd ${myCWD}/" >> mask_14_TKmerge_single_${RAW_DB%.db}.${slurmID}.plan
        echo "MARVEL TKmerge $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > mask_14_TKmerge_single_${RAW_DB%.db}.${slurmID}.version
        echo "DAZZLER Catrack $(git --git-dir=${DAZZLER_SOURCE_PATH}/DAZZ_DB/.git rev-parse --short HEAD)" >> mask_14_TKmerge_single_${RAW_DB%.db}.${slurmID}.version
    else 
        (>&2 echo "step ${currentStep} in RAW_REPMASK_TYPE ${RAW_REPMASK_TYPE} not supported")
        (>&2 echo "valid steps are: ${myTypes[${RAW_REPMASK_TYPE}]}")
        exit 1        
    fi    
else
    (>&2 echo "unknown RAW_REPMASK_TYPE ${RAW_REPMASK_TYPE}")
    (>&2 echo "supported types")
    x=0; while [ $x -lt ${#myTypes[*]} ]; do (>&2 echo "${myTypes[${x}]}"); done
    exit 1
fi

exit 0
