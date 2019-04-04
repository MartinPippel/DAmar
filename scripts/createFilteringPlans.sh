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

if [[ -z ${FIX_FILT_SCRUB_TYPE} ]]
then
	(>&2 echo "WARNING - Variable FIX_FILT_SCRUB_TYPE is not set. Use default mode: dalign!")
	FIX_FILT_SCRUB_TYPE=1
fi

if [[ -z "${PROJECT_ID}" ]]
then 
    (>&2 echo "ERROR - You have to specify a project id. Set variable PROJECT_ID")
    exit 1
fi

if [[ ! -n "${FIX_FILT_TYPE}" ]]
then 
    (>&2 echo "cannot create read patching scripts if variable FIX_FILT_TYPE is not set.")
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
    	
    if [[ -z ${FIX_FILT_OUTDIR} ]]
    then
        FIX_FILT_OUTDIR="m1"
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

    if [[ -n ${FIX_FILT_LAFILTER_DIF} && ${FIX_FILT_LAFILTER_DIF} -ne 0 ]]
    then
        FILT_LAFILTER_OPT="${FILT_LAFILTER_OPT} -d ${FIX_FILT_LAFILTER_DIF}"
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
    if [[ -n ${FIX_FILT_LAFILTER_MINTIPCOV} && ${FIX_FILT_LAFILTER_MINTIPCOV} -gt 0 ]]
    then
        FILT_LAFILTER_OPT="${FILT_LAFILTER_OPT} -z ${FIX_FILT_LAFILTER_MINTIPCOV}"
    fi            
    if [[ -n ${FIX_FILT_LAFILTER_MULTIMAPPER} && ${FIX_FILT_LAFILTER_MULTIMAPPER} -gt 0 ]]
    then
        FILT_LAFILTER_OPT="${FILT_LAFILTER_OPT} -w"
    fi
    if [[ -n ${FIX_FILT_LAFILTER_MAXREPEATMERGELEN} && ${FIX_FILT_LAFILTER_MAXREPEATMERGELEN} -gt 0 ]]
    then
        FILT_LAFILTER_OPT="${FILT_LAFILTER_OPT} -V ${FIX_FILT_LAFILTER_MAXREPEATMERGELEN}"
    fi
    if [[ -n ${FIX_FILT_LAFILTER_MAXREPEATMERGEWINDOW} && ${FIX_FILT_LAFILTER_MAXREPEATMERGEWINDOW} -gt 0 ]]
    then
    	FILT_LAFILTER_OPT="${FILT_LAFILTER_OPT} -W ${FIX_FILT_LAFILTER_MAXREPEATMERGEWINDOW}"
    fi
                
    if [[ -n ${FIX_FILT_LAFILTER_EXCLUDEREADS} || -n ${FIX_SCRUB_LAGAP_DISCARD_CHIMERS} ]]
    then
        if [[ -n ${FIX_SCRUB_LAGAP_DISCARD_CHIMERS} ]]
        then 
            ptype=""
            if [[ ${FIX_FILT_SCRUB_TYPE} -eq 1 ]]
            then 
                ptype="dalign"
            elif [[ ${FIX_FILT_SCRUB_TYPE} -eq 2 ]]
            then 
                ptype="repcomp"
            elif [[ ${FIX_FILT_SCRUB_TYPE} -eq 3 ]]
            then 
                ptype="forcealign"
            fi
            for x in $(seq 1 ${fixblocks})
            do     
                cat ${FIX_DB%.db}.${x}.${ptype}Gap.chimers.txt
            done > ${FIX_DB%.db}.${ptype}Gap.chimers.txt

            ## if additional reads have to be excluded append them to final gapChimer file
            if [[ -n ${FIX_FILT_LAFILTER_EXCLUDEREADS} ]]
            then 
                cat ${FIX_FILT_LAFILTER_EXCLUDEREADS} >> ${FIX_DB%.db}.${ptype}Gap.chimers.txt
            fi  
            FILT_LAFILTER_OPT="${FILT_LAFILTER_OPT} -x ${FIX_DB%.db}.${ptype}Gap.chimers.txt"
        else
            FILT_LAFILTER_OPT="${FILT_LAFILTER_OPT} -x ${FIX_FILT_LAFILTER_EXCLUDEREADS}"
        fi        
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
    if [[ -n ${FIX_FILT_LAFILTER_STITCH} && ${FIX_FILT_LAFILTER_STITCH} -gt 0 ]]
    then
        if [[ -n ${FIX_FILT_LAFILTER_STITCH_AGG} && ${FIX_FILT_LAFILTER_STITCH_AGG} -gt 0 ]]
        then
            FILT_LAFILTER_OPT="${FILT_LAFILTER_OPT} -S ${FIX_FILT_LAFILTER_STITCH}"
        else
            FILT_LAFILTER_OPT="${FILT_LAFILTER_OPT} -s ${FIX_FILT_LAFILTER_STITCH}"
        fi
    fi
    
    if [[ -n ${FIX_FILT_LAFILTER_TRIM} && ${FIX_FILT_LAFILTER_TRIM} -ne 0 ]]
    then
        if [[ -z ${SCRUB_LAQ_OPT} ]]
        then 
            setLAqOptions
        fi
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
    fi

    if [[ ${FIX_FILT_SCRUB_TYPE} -eq 1 ]]
    then 
        FIX_FILT_ENDING="dalignGap"
        if [[ -n ${FIX_FILT_LAFILTER_TRIM} && ${FIX_FILT_LAFILTER_TRIM} -ne 0 ]]
        then
            FILT_LAFILTER_OPT="${FILT_LAFILTER_OPT} -t trim1_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}_dalign -T"
        fi
        if [[ -n ${FIX_FILT_LAFILTER_REPEAT_IDX} ]]
        then 
        	if [[ ${#FIX_FILT_LAFILTER_REPEAT_IDX[@]} -eq 0 ]] ## its not a bash arroy
        	then
        		tmp=f$(echo ${SCRUB_LAREPEAT_OPT[${FIX_FILT_LAFILTER_REPEAT_IDX}]} | awk '{print $NF}')_dalign_${FIX_REPMASK_LAREPEAT_REPEATTRACK}
        	elif [[ ${#FIX_FILT_LAFILTER_REPEAT_IDX[@]} -eq 1 ]] ## bash arroy with one element
        	then
        		tmp=f$(echo ${SCRUB_LAREPEAT_OPT[${FIX_FILT_LAFILTER_REPEAT_IDX[0]}]} | awk '{print $NF}')_dalign_${FIX_REPMASK_LAREPEAT_REPEATTRACK}
        	elif [[ ${#FIX_FILT_LAFILTER_REPEAT_IDX[@]} -eq 2 ]] ## bash arroy with two element
        	then
        		
        		tmp1=f$(echo ${SCRUB_LAREPEAT_OPT[${FIX_FILT_LAFILTER_REPEAT_IDX[0]}]} | awk '{print $NF}')
        		tmp2=f$(echo ${SCRUB_LAREPEAT_OPT[${FIX_FILT_LAFILTER_REPEAT_IDX[1]}]} | awk '{print $NF}')
        		
        		if [[ ! -f .${FIX_DB%.db}.${tmp1}_${tmp2}_dalign_${FIX_REPMASK_LAREPEAT_REPEATTRACK}.d2 ]]
        		then
        			echo "WARNING missing track ${tmp1}_${tmp2}_dalign_${FIX_REPMASK_LAREPEAT_REPEATTRACK}!!!" 
        			echo "run TKcombine for tracks ${tmp1}_dalign_${FIX_REPMASK_LAREPEAT_REPEATTRACK} and ${tmp2}_dalign_${FIX_REPMASK_LAREPEAT_REPEATTRACK} on the fly"
        			
        			${MARVEL_PATH}/bin/TKcombine -v ${FIX_DB%.db} ${tmp1}_${tmp2}_dalign_${FIX_REPMASK_LAREPEAT_REPEATTRACK} ${tmp1}_dalign_${FIX_REPMASK_LAREPEAT_REPEATTRACK} ${tmp2}_dalign_${FIX_REPMASK_LAREPEAT_REPEATTRACK}
        			${MARVEL_PATH}/bin/TKcombine -v ${FIX_DB%.db} ${tmp1}_${tmp2}_dalign_${FIX_REPMASK_LAREPEAT_REPEATTRACK}_${FIX_REPMASK_TANMASK_TRACK}_dust ${FIX_REPMASK_TANMASK_TRACK} dust
        			
        			if [[ ! -f .${FIX_DB%.db}.${tmp1}_${tmp2}_dalign_${FIX_REPMASK_LAREPEAT_REPEATTRACK}.d2 ]]
        			then
        				(>&2 echo "ERROR could not create desired repeat track ${tmp1}_${tmp2}_dalign_${FIX_REPMASK_LAREPEAT_REPEATTRACK}!!!")
        				exit 1
        			fi  
        			      
        			tmp=${tmp1}_${tmp2}_dalign_${FIX_REPMASK_LAREPEAT_REPEATTRACK}			 	
        		fi 
        	else ## bash arroy with more than two element are not supported
        		(>&2 echo "More then two repeat tracks are nor supported yet!!")
        		exit 1
        	fi
            
            if [[ -n ${FIX_FILT_LAFILTER_DUST} ]]
            then 
            	FILT_LAFILTER_OPT="${FILT_LAFILTER_OPT} -D ${FIX_FILT_LAFILTER_DUST}"
            	FIX_FILT_LAFILTER_REPEATTRACK=${tmp}
        	else
        		FIX_FILT_LAFILTER_REPEATTRACK=${tmp}_${FIX_REPMASK_TANMASK_TRACK}_dust
			fi
            FILT_LAFILTER_OPT="${FILT_LAFILTER_OPT} -r ${FIX_FILT_LAFILTER_REPEATTRACK}"
        fi
    elif [[ ${FIX_FILT_SCRUB_TYPE} -eq 2 ]]
    then 
        FIX_FILT_ENDING="repcompGap"  
        if [[ -n ${FIX_FILT_LAFILTER_TRIM} && ${FIX_FILT_LAFILTER_TRIM} -ne 0 ]]
        then
            FILT_LAFILTER_OPT="${FILT_LAFILTER_OPT} -t trim1_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}_repcomp -T"
        fi          
        if [[ -n ${FIX_FILT_LAFILTER_REPEAT_IDX} ]]
        then 
        	if [[ ${#FIX_FILT_LAFILTER_REPEAT_IDX[@]} -eq 0 ]] ## its not a bash arroy
        	then
        		tmp=f$(echo ${SCRUB_LAREPEAT_OPT[${FIX_FILT_LAFILTER_REPEAT_IDX}]} | awk '{print $NF}')_repcomp_${FIX_REPMASK_LAREPEAT_REPEATTRACK}
        	elif [[ ${#FIX_FILT_LAFILTER_REPEAT_IDX[@]} -eq 1 ]] ## bash arroy with one element
        	then
        		tmp=f$(echo ${SCRUB_LAREPEAT_OPT[${FIX_FILT_LAFILTER_REPEAT_IDX[0]}]} | awk '{print $NF}')_repcomp_${FIX_REPMASK_LAREPEAT_REPEATTRACK}
        	elif [[ ${#FIX_FILT_LAFILTER_REPEAT_IDX[@]} -eq 2 ]] ## bash arroy with two element
        	then
        		
        		tmp1=f$(echo ${SCRUB_LAREPEAT_OPT[${FIX_FILT_LAFILTER_REPEAT_IDX[0]}]} | awk '{print $NF}')
        		tmp2=f$(echo ${SCRUB_LAREPEAT_OPT[${FIX_FILT_LAFILTER_REPEAT_IDX[1]}]} | awk '{print $NF}')
        		
        		if [[ ! -f .${FIX_DB%.db}.${tmp1}_${tmp2}_repcomp_${FIX_REPMASK_LAREPEAT_REPEATTRACK}.d2 ]]
        		then
        			echo "WARNING missing track ${tmp1}_${tmp2}_repcomp_${FIX_REPMASK_LAREPEAT_REPEATTRACK}!!!" 
        			echo "run TKcombine for tracks ${tmp1}_repcomp_${FIX_REPMASK_LAREPEAT_REPEATTRACK} and ${tmp2}_repcomp_${FIX_REPMASK_LAREPEAT_REPEATTRACK} on the fly"
        			
        			${MARVEL_PATH}/bin/TKcombine -v ${FIX_DB%.db} ${tmp1}_${tmp2}_repcomp_${FIX_REPMASK_LAREPEAT_REPEATTRACK} ${tmp1}_repcomp_${FIX_REPMASK_LAREPEAT_REPEATTRACK} ${tmp2}_repcomp_${FIX_REPMASK_LAREPEAT_REPEATTRACK}
        			${MARVEL_PATH}/bin/TKcombine -v ${FIX_DB%.db} ${tmp1}_${tmp2}_repcomp_${FIX_REPMASK_LAREPEAT_REPEATTRACK}_${FIX_REPMASK_TANMASK_TRACK}_dust ${FIX_REPMASK_TANMASK_TRACK} dust
        			
        			if [[ ! -f .${FIX_DB%.db}.${tmp1}_${tmp2}_repcomp_${FIX_REPMASK_LAREPEAT_REPEATTRACK}.d2 ]]
        			then
        				(>&2 echo "ERROR could not create desired repeat track ${tmp1}_${tmp2}_repcomp_${FIX_REPMASK_LAREPEAT_REPEATTRACK}!!!")
        				exit 1
        			fi  
        			      
        			tmp=${tmp1}_${tmp2}_repcomp_${FIX_REPMASK_LAREPEAT_REPEATTRACK}			 	
        		fi 
        	else ## bash arroy with more than two element are not supported
        		(>&2 echo "More then two repeat tracks are nor supported yet!!")
        		exit 1
        	fi
            if [[ -n ${FIX_FILT_LAFILTER_DUST} ]]
            then 
            	FILT_LAFILTER_OPT="${FILT_LAFILTER_OPT} -D ${FIX_FILT_LAFILTER_DUST}"
            	FIX_FILT_LAFILTER_REPEATTRACK=${tmp}
        	else
        		FIX_FILT_LAFILTER_REPEATTRACK=${tmp}_${FIX_REPMASK_TANMASK_TRACK}_dust
			fi            
            FILT_LAFILTER_OPT="${FILT_LAFILTER_OPT} -r ${FIX_FILT_LAFILTER_REPEATTRACK}"
        fi
    elif [[ ${FIX_FILT_SCRUB_TYPE} -eq 3 ]]
    then 
        FIX_FILT_ENDING="forcealignGap"
        if [[ -n ${FIX_FILT_LAFILTER_TRIM} && ${FIX_FILT_LAFILTER_TRIM} -ne 0 ]]
        then
            FILT_LAFILTER_OPT="${FILT_LAFILTER_OPT} -t trim1_d${FIX_SCRUB_LAQ_QTRIMCUTOFF}_s${FIX_SCRUB_LAQ_MINSEG}_forcealign -T"
        fi
        if [[ -n ${FIX_FILT_LAFILTER_REPEAT_IDX} ]]
        then 
        	if [[ ${#FIX_FILT_LAFILTER_REPEAT_IDX[@]} -eq 0 ]] ## its not a bash arroy
        	then
        		tmp=f$(echo ${SCRUB_LAREPEAT_OPT[${FIX_FILT_LAFILTER_REPEAT_IDX}]} | awk '{print $NF}')_forcealign_${FIX_REPMASK_LAREPEAT_REPEATTRACK}
        	elif [[ ${#FIX_FILT_LAFILTER_REPEAT_IDX[@]} -eq 1 ]] ## bash arroy with one element
        	then
        		tmp=f$(echo ${SCRUB_LAREPEAT_OPT[${FIX_FILT_LAFILTER_REPEAT_IDX[0]}]} | awk '{print $NF}')_forcealign_${FIX_REPMASK_LAREPEAT_REPEATTRACK}
        	elif [[ ${#FIX_FILT_LAFILTER_REPEAT_IDX[@]} -eq 2 ]] ## bash arroy with two element
        	then
        		
        		tmp1=f$(echo ${SCRUB_LAREPEAT_OPT[${FIX_FILT_LAFILTER_REPEAT_IDX[0]}]} | awk '{print $NF}')
        		tmp2=f$(echo ${SCRUB_LAREPEAT_OPT[${FIX_FILT_LAFILTER_REPEAT_IDX[1]}]} | awk '{print $NF}')
        		
        		if [[ ! -f .${FIX_DB%.db}.${tmp1}_${tmp2}_forcealign_${FIX_REPMASK_LAREPEAT_REPEATTRACK}.d2 ]]
        		then
        			echo "WARNING missing track ${tmp1}_${tmp2}_forcealign_${FIX_REPMASK_LAREPEAT_REPEATTRACK}!!!" 
        			echo "run TKcombine for tracks ${tmp1}_forcealign_${FIX_REPMASK_LAREPEAT_REPEATTRACK} and ${tmp2}_forcealign_${FIX_REPMASK_LAREPEAT_REPEATTRACK} on the fly"
        			
        			${MARVEL_PATH}/bin/TKcombine -v ${FIX_DB%.db} ${tmp1}_${tmp2}_forcealign_${FIX_REPMASK_LAREPEAT_REPEATTRACK} ${tmp1}_forcealign_${FIX_REPMASK_LAREPEAT_REPEATTRACK} ${tmp2}_forcealign_${FIX_REPMASK_LAREPEAT_REPEATTRACK}
        			${MARVEL_PATH}/bin/TKcombine -v ${FIX_DB%.db} ${tmp1}_${tmp2}_forcealign_${FIX_REPMASK_LAREPEAT_REPEATTRACK}_${FIX_REPMASK_TANMASK_TRACK}_dust ${FIX_REPMASK_TANMASK_TRACK} dust
        			
        			if [[ ! -f .${FIX_DB%.db}.${tmp1}_${tmp2}_forcealign_${FIX_REPMASK_LAREPEAT_REPEATTRACK}.d2 ]]
        			then
        				(>&2 echo "ERROR could not create desired repeat track ${tmp1}_${tmp2}_forcealign_${FIX_REPMASK_LAREPEAT_REPEATTRACK}!!!")
        				exit 1
        			fi  
        			      
        			tmp=${tmp1}_${tmp2}_forcealign_${FIX_REPMASK_LAREPEAT_REPEATTRACK}			 	
        		fi 
        	else ## bash arroy with more than two element are not supported
        		(>&2 echo "More then two repeat tracks are nor supported yet!!")
        		exit 1
        	fi
            if [[ -n ${FIX_FILT_LAFILTER_DUST} ]]
            then 
            	FILT_LAFILTER_OPT="${FILT_LAFILTER_OPT} -D ${FIX_FILT_LAFILTER_DUST}"
            	FIX_FILT_LAFILTER_REPEATTRACK=${tmp}
        	else
        		FIX_FILT_LAFILTER_REPEATTRACK=${tmp}_${FIX_REPMASK_TANMASK_TRACK}_dust
			fi

            FILT_LAFILTER_OPT="${FILT_LAFILTER_OPT} -r ${FIX_FILT_LAFILTER_REPEATTRACK}"
        fi
    else
        (>&2 echo "step ${currentStep} in FIX_FILT_SCRUB_TYPE ${FIX_FILT_SCRUB_TYPE} not supported")
        exit 1
	fi
}

function setTKmergeOptions() 
{
    FILT_TKMERGE_OPT=""
    if [[ -n ${FIX_FILT_TKMERGE_DELETE} && ${FIX_FILT_TKMERGE_DELETE} -ne 0 ]]
    then
        FILT_TKMERGE_OPT="${FILT_TKMERGE_OPT} -d"
    fi
}

function setLAmergeOptions()
{
    FILT_LAMERGE_OPT=""
    if [[ -n ${FIX_FILT_LAMERGE_NFILES} && ${FIX_FILT_LAMERGE_NFILES} -gt 0 ]]
    then
        FILT_LAMERGE_OPT="${FILT_LAMERGE_OPT} -n ${FIX_FILT_LAMERGE_NFILES}"
    fi
}

## ensure some paths
if [[ -z "${MARVEL_SOURCE_PATH}" || ! -d "${MARVEL_SOURCE_PATH}" ]]
then 
    (>&2 echo "ERROR - You have to set MARVEL_SOURCE_PATH. Used to report git version.")
    exit 1
fi

myTypes=("1-createSubdirFILT_FSUFFIX, 2-LAfilter, 3-LAmerge")
#type-0 steps: 1-createSubdirFILT_FSUFFIX, 2-LAfilter, 3-LAmerge
if [[ ${FIX_FILT_TYPE} -eq 0 ]]
then 
    ### create sub-directory and link relevant DB and Track files
    if [[ ${currentStep} -eq 1 ]]
    then
        ### clean up plans 
        for x in $(ls filt_01_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 

        setLAfilterOptions        

        echo "if [[ -d ${FIX_FILT_OUTDIR} ]]; then mv ${FIX_FILT_OUTDIR} ${FIX_FILT_OUTDIR}_$(date '+%Y-%m-%d_%H-%M-%S'); fi && mkdir ${FIX_FILT_OUTDIR} && ln -s -r .${FIX_DB%db}.* ${FIX_DB%db}.db ${FIX_FILT_OUTDIR}" > filt_01_createSubDir_single_${FIX_DB%.db}.${slurmID}.plan
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > filt_01_createSubDir_single_${FIX_DB%.db}.${slurmID}.version         
    ### LAfilter
    elif [[ ${currentStep} -eq 2 ]]
    then
        ### clean up plans 
        for x in $(ls filt_02_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
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
                if [[ -f filt.round${rnd}_02_LAfilter_block_${FIX_DB%.db}.${slurmID}.plan ]]
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
                    
                    echo "${MARVEL_PATH}/bin/LAfilter${FILT_LAFILTER_OPT}${addOpt} ${FIX_FILT_OUTDIR}/${FIX_DB%.db} ${FIX_DB%.db}.${x}.${FIX_FILT_ENDING}.las ${FIX_FILT_OUTDIR}/${FIX_DB%.db}.${x}.filt_R1.las"
                done > filt.round1_02_LAfilter_block_${FIX_DB%.db}.${slurmID}.plan 
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
                done > filt_02_LAfilter_block_${FIX_DB%.db}.${slurmID}.plan 
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
                done > filt.round$((${rnd}+1))_02_LAfilter_block_${FIX_DB%.db}.${slurmID}.plan 
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
                echo "${MARVEL_PATH}/bin/LAfilter${FILT_LAFILTER_OPT}${addOpt} ${FIX_FILT_OUTDIR}/${FIX_DB%.db} ${FIX_DB%.db}.${x}.${FIX_FILT_ENDING}.las ${FIX_FILT_OUTDIR}/${FIX_DB%.db}.${x}.filt.las"
            done > filt_02_LAfilter_block_${FIX_DB%.db}.${slurmID}.plan 
        fi    
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > filt_02_LAfilter_block_${FIX_DB%.db}.${slurmID}.version
    #### LAmerge
    elif [[ ${currentStep} -eq 3 ]]
    then
        ### clean up plans 
        for x in $(ls filt_03_*_*_${FIX_DB%.db}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done 
        
        if [[ -z ${FILT_LAFILTER_OPT} ]]
        then 
            setLAfilterOptions
        fi
        ### find and set LAmerge options 
        setLAmergeOptions
        
        echo "${MARVEL_PATH}/bin/LAmerge${FILT_LAMERGE_OPT} -S filt ${FIX_FILT_OUTDIR}/${FIX_DB%.db} ${FIX_FILT_OUTDIR}/${FIX_DB%.db}.filt.las" > filt_03_LAmerge_single_${FIX_DB%.db}.${slurmID}.plan
        echo "MARVEL $(git --git-dir=${MARVEL_SOURCE_PATH}/.git rev-parse --short HEAD)" > filt_03_LAmerge_single_${FIX_DB%.db}.${slurmID}.version             
    else
        (>&2 echo "step ${currentStep} in FIX_FILT_TYPE ${FIX_FILT_TYPE} not supported")
        (>&2 echo "valid steps are: ${myTypes[${FIX_FILT_TYPE}]}")
        exit 1            
    fi    
else
    (>&2 echo "unknown FIX_FILT_TYPE ${FIX_FILT_TYPE}")
    (>&2 echo "supported types")
    x=0; while [ $x -lt ${#myTypes[*]} ]; do (>&2 echo "${myTypes[${x}]}"); done
    exit 1
fi

exit 0
