#!/bin/bash 

configFile=$1
currentPhase=$2
currentStep=$3
slurmID=$4
if [[ -n $5 ]]
then 
	resumeIdx=$5
else
	resumeIdx=0
fi

retrySubmit=3

echo "createAndSubmitMarvelSlurmJobs.sh ${configFile} ${currentPhase} ${currentStep} ${slurmID}"
function isNumber 
{
  re='^[0-9]+$'
  if ! [[ $1 =~ $re ]]
  then
   exit 1
  fi
}

function prependLeadingNull
{   
 if [[ ! -n ${1} ]]
 then
   echo "prependLeadingNull needs an argument!!!" 
   exit 1
 fi

 if [[ ${1} -lt 10 ]]
 then
  echo "0${1}"
 else
  echo "$1"
 fi
}

source ${configFile}

function getPhaseFilePrefix()
{
	if [[ ${currentPhase} -eq -2 ]]
    then
        echo "qc"
	elif [[ ${currentPhase} -eq -1 ]]
    then
        echo "mito"
	elif [[ ${currentPhase} -eq 0 ]]
    then
        echo "cover"	
    elif [[ ${currentPhase} -eq 1 || ${currentPhase} -eq 3 ]]
    then
        echo "mask"
    elif [[ ${currentPhase} -eq 2 ]]
    then
        echo "fix"
    elif [[ ${currentPhase} -eq 4 ]]
    then
        echo "scrub"
    elif [[ ${currentPhase} -eq 5 ]]
    then
        if [[ ${FIX_FILT_TYPE} -eq 0 && ${currentStep} -eq 2 && -n ${FIX_FILT_LAFILTER_RMSYMROUNDS} && ${FIX_FILT_LAFILTER_RMSYMROUNDS} -gt 0 ]]
        then
            if [[ -f filt_02_LAfilter_block_${FIX_DB%.db}.${slurmID}.plan ]]
            then
                 echo "filt"
            else        
                ## check what is the current round
                for rnd in $(seq ${FIX_FILT_LAFILTER_RMSYMROUNDS} -1 0)
                do
                    if [[ -f filt.round${rnd}_02_LAfilter_block_${FIX_DB%.db}.${slurmID}.plan ]]
                    then
                        break;
                    fi
                done
                echo "filt.round${rnd}"
            fi
        else     
            echo "filt"
        fi
    elif [[ ${currentPhase} -eq 6 ]]
    then
        echo "tour"
    elif [[ ${currentPhase} -eq 7 ]]
    then
        echo "corr"
    elif [[ ${currentPhase} -eq 8 ]]
    then
        echo "cont"
    elif [[ ${currentPhase} -eq 9 ]]
    then
        echo "arrow"
    elif [[ ${currentPhase} -eq 10 ]]
    then
        echo "purgeHaplotigs"       
    elif [[ ${currentPhase} -eq 11 ]]
    then
        echo "freebayes"                     
    elif [[ ${currentPhase} -eq 12 ]]
    then
        echo "phase"                     
    elif [[ ${currentPhase} -eq 13 ]]
    then
        echo "10x"                     
    elif [[ ${currentPhase} -eq 14 ]]
    then
        echo "bionano"
    elif [[ ${currentPhase} -eq 15 ]]
    then
        echo "hic"                     
    else
        (>&2 echo "unknown MARVEL phase ${currentPhase}! Supported values (-2-qc, -1-mito, 1-mask, 2-fix, 3-mask, 4-scrub, 5-filt, 6-tour, 7-corr, 8-cont, 9-arrow, 10-purgeHaplotigs, 11-freebayes, 12-phase, 13-10x, 14-bionano, 15-hic)")
        exit 1
    fi
}

function getCurrentDB()
{	
	if [[ ${currentPhase} -lt 0 ]]
    then 
        echo "${RAW_DB%.db}"
	elif [[ ${currentPhase} -eq 0 ]]
    then 
        echo "${RAW_DAZZ_DB%.db}"
    elif [[ ${currentPhase} -lt 3 ]]
    then 
        echo "${RAW_DB%.db}"
	elif [[ ${currentPhase} -lt 8 ]]
    then 
       	echo "${FIX_DB%.db}"       
    else
        echo "${CONT_DB%.db}"
    fi
}
echo $(pwd)

if [[ ${resumeIdx} -eq 0 ]]
then
	### create current plan 
	${SUBMIT_SCRIPTS_PATH}/createCommandPlan.sh ${configFile} ${currentPhase} ${currentStep} ${slurmID}
	if [ $? -ne 0 ]
	then 
    	(>&2 echo "createCommandPlan.sh failed some how. Stop here.")
    	exit 1      
	fi 
fi
 
db=$(getCurrentDB)
prefix=$(getPhaseFilePrefix)
cjobid=$(prependLeadingNull ${currentStep})
if ! ls ${prefix}_${cjobid}_*_*_${db}.${slurmID}.plan 1> /dev/null 2>&1;
then
    (>&2 echo "createAndSubmitMarvelSlumJobs.sh - missing file: ${prefix}_${cjobid}_*_*_${db}.${slurmID}.plan")
    exit 1
fi
### get job name 
jname=$(ls ${prefix}_${cjobid}_*_*_${db}.${slurmID}.plan | awk -F \_ '{print $3}')
jtype=$(ls ${prefix}_${cjobid}_*_*_${db}.${slurmID}.plan | awk -F \_ '{print $4}')
if [[ ${jtype} != "block" && ${jtype} != "single" ]]
then
    (>&2 echo "unknown slurm type: ${jtype}. valid types: block, single")
    exit 1
fi

### create submit scripts

if [[ ${resumeIdx} -eq 0 ]]
then
	### setup runtime conditions, time, memory, etc 
	MEM=${MEM_DEFAULT}
	TMP="MEM_${jname}"
	if [[ -n ${!TMP} ]]
	then 
	    MEM=${!TMP}
	fi
	TIME=${TIME_DEFAULT}
	TMP="TIME_${jname}"
	if [[ -n ${!TMP} ]]
	then
	    TIME=${!TMP}
	fi
	CORES=${THREADS_DEFAULT}
	TMP="THREADS_${jname}"
	if [[ -n ${!TMP} ]]
	then
	    CORES=${!TMP}
	fi
	NTASKS_PER_NODE=""
	TMP="TASKS_${jname}"
	if [[ -n ${!TMP} ]]
	then
	    NTASKS_PER_NODE="${!TMP}"
	fi
	STEPSIZE=""
	TMP="STEPSIZE_${jname}"
	if [[ -n ${!TMP} ]]
	then
	    STEPSIZE=":${!TMP}"
	fi
	TMP="PARTITION_${jname}"
	if [[ -n ${!TMP} ]]
	then
	    SLURM_PARTITION="${!TMP}"
	fi
	
	JOBS=$(wc -l ${prefix}_${cjobid}_${jname}_${jtype}_${db}.${slurmID}.plan | awk '{print $1}')
	log_folder=log_${prefix}_${jname}_${db}
	mkdir -p ${log_folder}
	first=1
	if [[ ${JOBS} -gt 9999 ]]
	then
	    from=1
	    to=9999
	    d=1
	    while [[ $from -lt ${JOBS} ]]
	    do
	        file=${prefix}_${cjobid}_${jname}_${jtype}_${db}.${slurmID}.${d}
	        sed -n ${from},${to}p ${prefix}_${cjobid}_${jname}_${jtype}_${db}.${slurmID}.plan > ${file}.plan
	        jobs=$((${to}-${from}+1))
	        ### create slurm submit file
	        if [[ ${jtype} == "block" ]]
	        then 
	            echo "#!/bin/bash
#SBATCH -J ${PROJECT_ID}_p${currentPhase}s${currentStep}
#SBATCH -p ${SLURM_PARTITION}
#SBATCH -a 1-${jobs}${STEPSIZE}
#SBATCH -c ${CORES} # Number of cores 
#SBATCH -n 1 # number of nodes
#SBATCH -o ${log_folder}/${prefix}_${cjobid}_${db}_${d}_%A_%a.out # Standard output 
#SBATCH -e ${log_folder}/${prefix}_${cjobid}_${db}_${d}_%A_%a.err # Standard error
#SBATCH --time=${TIME}
#SBATCH --mem-per-cpu=$((${MEM}/${CORES}))
#SBATCH --mail-user=pippel@mpi-cbg.de
#SBATCH --mail-type=FAIL" > ${file}.slurm
	            if [[ -n ${NTASKS_PER_NODE} ]]
	            then
	                echo "#SBATCH --ntasks-per-node=${NTASKS_PER_NODE}" >> ${file}.slurm
	            fi 
	        	if [[ -n ${SLURM_NUMACTL} && ${SLURM_NUMACTL} -gt 0  ]]
				then	
					echo -e "#SBATCH --mem_bind=verbose,local" >> ${file}.slurm			
				fi
                if [[ -n ${SLURM_ACCOUNT} ]]
                then
                    echo "#SBATCH -A ${SLURM_ACCOUNT}" >> ${file}.slurm
                fi
				if [[ ${prefix} == "arrow" || ${prefix} == "freebayes" || ${prefix} == "hic" || ${prefix} == "phase" || ${prefix} == "qc" ]]
				then
					echo -e "\n${PACBIO_BASE_ENV}" >> ${file}.slurm
				elif [[ ${prefix} == "purgeHaplotigs" ]]
				then	
					echo -e "\n${PURGEHAPLOTIGS_ENV}" >> ${file}.slurm
				fi
				
	            echo "export PATH=${MARVEL_PATH}/bin:\$PATH
export PATH=${MARVEL_PATH}/scripts:\$PATH
export PYTHONPATH=${MARVEL_PATH}/lib.python:\$PYTHONPATH

FIRSTJOB=0
LASTJOB=\$(wc -l ${file}.plan | awk '{print \$1}')

beg=\$(date +%s)
echo \"${file}.plan beg $beg\"

i=0;
while [[ \$i -lt \$SLURM_ARRAY_TASK_STEP ]]
do
  index=\$((\$SLURM_ARRAY_TASK_ID+\$i+\${FIRSTJOB}))
  echo \"i \$i index: \$index\"
  if [[ \$index -le \$LASTJOB ]]
  then
    echo \"eval line \$index\"
    eval \$(sed -n \${index}p ${file}.plan) || exit 100
  fi
  i=\$((\$i+1))
done

end=\$(date +%s)
echo \"${file}.plan end \$end\"
echo \"${file}.plan run time: \$((\${end}-\${beg}))\"" >> ${file}.slurm
			else    # jtype == single ? can this really happen > 9999 jobs that are sequentially executed 
	            echo "#!/bin/bash
#SBATCH -J ${PROJECT_ID}_p${currentPhase}s${currentStep}
#SBATCH -p ${SLURM_PARTITION}
#SBATCH -c ${CORES} # Number of cores
#SBATCH -n 1 # number of nodes
#SBATCH -o ${log_folder}/${prefix}_${cjobid}_${db}_${d}_%A.out # Standard output
#SBATCH -e ${log_folder}/${prefix}_${cjobid}_${db}_${d}_%A.err # Standard error
#SBATCH --time=${TIME}
#SBATCH --mem-per-cpu=$((${MEM}/${CORES}))
#SBATCH --mail-user=pippel@mpi-cbg.de
#SBATCH --mail-type=FAIL" > ${file}.slurm
			if [[ -n ${SLURM_NUMACTL} && ${SLURM_NUMACTL} -gt 0  ]]
			then	
				echo -e "#SBATCH --mem_bind=verbose,local" >> ${file}.slurm			
			fi
		    if [[ -n ${SLURM_ACCOUNT} ]]
            then
                echo "#SBATCH -A ${SLURM_ACCOUNT}" >> ${file}.slurm
            fi

			if [[ ${prefix} == "arrow" || ${prefix} == "freebayes" || ${prefix} == "hic" || ${prefix} == "phase" || ${prefix} == "qc" ]]
			then
				echo -e "\n${PACBIO_BASE_ENV}" >> ${file}.slurm
			elif [[ ${prefix} == "purgeHaplotigs" ]]
			then	
				echo -e "\n${PURGEHAPLOTIGS_ENV}" >> ${file}.slurm
			fi

			echo "export PATH=${MARVEL_PATH}/bin:\$PATH
export PATH=${MARVEL_PATH}/scripts:\$PATH
export PYTHONPATH=${MARVEL_PATH}/lib.python:\$PYTHONPATH

FIRSTJOB=1
LASTJOB=\$(wc -l ${file}.plan | awk '{print \$1}')

beg=\$(date +%s)
echo \"${file}.plan beg \$beg\"

i=\${FIRSTJOB}
while [[ \$i -le \${LASTJOB} ]]
do 
  echo \"eval line \$i\"
  eval \$(sed -n \${i}p ${file}.plan) || exit 100
  i=\$((\$i+1))
done

end=\$(date +%s)
echo \"${file}.plan end \$end\"
echo \"${file}.plan run time: $((${end}-${beg}))\"" >> ${file}}.slurm
	        fi
	        d=$(($d+1))
	        from=$((${to}+1))
	        to=$((${to}+9999))
	        if [[ $to -gt ${JOBS} ]]
	        then
	            to=${JOBS}
	        fi
	        sleep 5
	    done
	else ## less then 9999 jobs 
	    jobs=${JOBS}
	    file=${prefix}_${cjobid}_${jname}_${jtype}_${db}.${slurmID}
	    ### create slurm submit file
	    if [[ ${jtype} == "block" ]]
	    then 
	        echo "#!/bin/bash
#SBATCH -J ${PROJECT_ID}_p${currentPhase}s${currentStep}
#SBATCH -p ${SLURM_PARTITION}
#SBATCH -a 1-${jobs}${STEPSIZE}
#SBATCH -c ${CORES} # Number of cores 
#SBATCH -n 1 # number of nodes
#SBATCH -o ${log_folder}/${prefix}_${cjobid}_${db}_%A_%a.out # Standard output
#SBATCH -e ${log_folder}/${prefix}_${cjobid}_${db}_%A_%a.err # Standard error
#SBATCH --time=${TIME}
#SBATCH --mem-per-cpu=$((${MEM}/${CORES}))
#SBATCH --mail-user=pippel@mpi-cbg.de
#SBATCH --mail-type=FAIL" > ${file}.slurm
	        if [[ -n ${NTASKS_PER_NODE} ]]
	        then
	            echo "#SBATCH --ntasks-per-node=${NTASKS_PER_NODE}" >> ${file}.slurm
	        fi 
	        if [[ -n ${SLURM_NUMACTL} && ${SLURM_NUMACTL} -gt 0  ]]
			then	
				echo -e "#SBATCH --mem_bind=verbose,local" >> ${file}.slurm			
			fi
            if [[ -n ${SLURM_ACCOUNT} ]]
            then
                echo "#SBATCH -A ${SLURM_ACCOUNT}" >> ${file}.slurm
            fi	        
	        
			if [[ ${prefix} == "arrow" || ${prefix} == "freebayes" || ${prefix} == "hic" || ${prefix} == "phase" || ${prefix} == "qc" ]]
			then
				echo -e "\n${PACBIO_BASE_ENV}" >> ${file}.slurm
			elif [[ ${prefix} == "purgeHaplotigs" ]]
			then	
				echo -e "\n${PURGEHAPLOTIGS_ENV}" >> ${file}.slurm
			fi
	
	        echo "export PATH=${MARVEL_PATH}/bin:\$PATH
export PATH=${MARVEL_PATH}/scripts:\$PATH
export PYTHONPATH=${MARVEL_PATH}/lib.python:\$PYTHONPATH

FIRSTJOB=0
LASTJOB=\$(wc -l ${file}.plan | awk '{print \$1}')

beg=\$(date +%s)
echo \"${file}.plan beg $beg\"

i=0;
while [ \$i -lt \$SLURM_ARRAY_TASK_STEP ]
do
  index=\$((\$SLURM_ARRAY_TASK_ID+\$i+\${FIRSTJOB}))
  echo \"i \$i index: \$index\"
  if [[ \$index -le \$LASTJOB ]]
  then
    echo \"eval line \$index\"
    eval \$(sed -n \${index}p ${file}.plan) || exit 100
  fi
  i=\$((\$i+1))
done

end=\$(date +%s)
echo \"${file}.plan end \$end\"
echo \"${file}.plan run time: \$((\${end}-\${beg}))\"" >> ${file}.slurm
	    else
	        echo "#!/bin/bash
#SBATCH -J ${PROJECT_ID}_p${currentPhase}s${currentStep}
#SBATCH -p ${SLURM_PARTITION}
#SBATCH -c ${CORES} # Number of cores
#SBATCH -n 1 # number of nodes
#SBATCH -o ${log_folder}/${prefix}_${cjobid}_${db}_%A.out # Standard output
#SBATCH -e ${log_folder}/${prefix}_${cjobid}_${db}_%A.err # Standard error
#SBATCH --time=${TIME}
#SBATCH --mem-per-cpu=$((${MEM}/${CORES}))
#SBATCH --mail-user=pippel@mpi-cbg.de
#SBATCH --mail-type=FAIL" > ${file}.slurm

			if [[ -n ${SLURM_NUMACTL} && ${SLURM_NUMACTL} -gt 0  ]]
			then	
				echo -e "#SBATCH --mem_bind=verbose,local" >> ${file}.slurm			
			fi
            if [[ -n ${SLURM_ACCOUNT} ]]
            then
                echo "#SBATCH -A ${SLURM_ACCOUNT}" >> ${file}.slurm
            fi
			if [[ ${prefix} == "arrow" || ${prefix} == "freebayes" || ${prefix} == "hic" || ${prefix} == "phase" || ${prefix} == "qc" ]]
			then
				echo -e "\n${PACBIO_BASE_ENV}" >> ${file}.slurm
			elif [[ ${prefix} == "purgeHaplotigs" ]]
			then	
				echo -e "\n${PURGEHAPLOTIGS_ENV}" >> ${file}.slurm
			fi

			echo "export PATH=${MARVEL_PATH}/bin:\$PATH
export PATH=${MARVEL_PATH}/scripts:\$PATH
export PYTHONPATH=${MARVEL_PATH}/lib.python:\$PYTHONPATH

FIRSTJOB=1
LASTJOB=\$(wc -l ${file}.plan | awk '{print \$1}')

beg=\$(date +%s)
echo \"${file}.plan beg \$beg\"

i=\${FIRSTJOB};
while [[ \$i -le \${LASTJOB} ]]
do
  echo \"eval line \$i\"
  eval \$(sed -n \${i}p ${file}.plan) || exit 100
  i=\$((\$i+1))
done

end=\$(date +%s)
echo \"${file}.plan end \$end\"
echo \"${file}.plan run time: \$((\${end}-\${beg}))\"" >> ${file}.slurm
	    fi
	fi
fi

if [[ ${resumeIdx} -eq 0 ]]
then
	if [[ ${JOBS} -gt 9999 ]]
	then
		resumeIdx=1
		file=${prefix}_${cjobid}_${jname}_${jtype}_${db}.${slurmID}.${resumeIdx}
	else
		file=${prefix}_${cjobid}_${jname}_${jtype}_${db}.${slurmID}
	fi		
else
	file=${prefix}_${cjobid}_${jname}_${jtype}_${db}.${slurmID}.${resumeIdx}
fi

retry=0
TMPRET=-1
wait=120
while [[ "${TMPRET}" == "-1" && ${retry} -lt ${retrySubmit} ]]
do
	if [[ ${retry} -gt 0 ]]
	then
		echo "try to restart job ${file}.slurm ${retry}/${retrySubmit} - wait $((${retry}*${wait})) seconds"
		sleep $((${retry}*${wait}))
	fi
	TMPRET=$(sbatch ${file}.slurm) && isNumber ${TMPRET##* } || TMPRET=-1            		
	retry=$((${retry}+1))
done

if [[ "${TMPRET}" == "-1" ]]
then
	(>&2 echo "Unable to submit job ${file}.slurm. Stop here.")
	exit 1
fi
echo "submit ${file}.slurm ${TMPRET##* }"
RET="${TMPRET##* }"

 
foundNext=0 
### add if account is necessary
appAccount=""
if [[ -n ${SLURM_ACCOUNT} ]]
then
	appAccount=" -A ${SLURM_ACCOUNT}"
fi

if [[ ${resumeIdx} -gt 0 ]]
then 
	if [[ -f ${prefix}_${cjobid}_${jname}_${jtype}_${db}.${slurmID}.$((${resumeIdx}+1)).slurm ]]
	then 
		sbatch${appAccount} --job-name=${PROJECT_ID}_p${currentPhase}s${currentStep+1} -o ${prefix}_step${currentStep}_${RAW_DB%.db}_${slurmID}.out -e ${prefix}_step${currentStep}_${RAW_DB%.db}_${slurmID}.err -n1 -c1 -p ${SLURM_PARTITION} --time=01:00:00 --mem-per-cpu=6g --dependency=afterok:${RET##* } --wrap="bash ${SUBMIT_SCRIPTS_PATH}/createAndSubmitMarvelSlurmJobs.sh ${configFile} ${currentPhase} ${currentStep} $slurmID $((${resumeIdx}+1))"
		foundNext=1
	fi	
fi

if [[ ${currentPhase} -eq -2 ]]
then
	if [[ $((${currentStep}+1)) -le ${RAW_QC_SUBMIT_SCRIPTS_TO} ]]
    then
    	sbatch${appAccount} --job-name=${PROJECT_ID}_p${currentPhase}s$((${currentStep+1})) -o ${prefix}_step$((${currentStep}+1))_${db}_${slurmID}.out -e ${prefix}_step$((${currentStep}+1))_${db%.db}_${slurmID}.err -n1 -c1 -p ${SLURM_PARTITION} --time=01:00:00 --mem-per-cpu=6g --dependency=afterok:${RET##* } --wrap="bash ${SUBMIT_SCRIPTS_PATH}/createAndSubmitMarvelSlurmJobs.sh ${configFile} ${currentPhase} $((${currentStep}+1)) $slurmID"
        foundNext=1
	else
		### we have to create the coverage directory and change into that dir 
		if [[ $(echo "$(pwd)" | awk -F \/ '{print $NF}') -eq ${MASH_DIR} ]]
		then
			if [[ ${RAW_MITO_SUBMIT_SCRIPTS_FROM} -gt 0 && ${RAW_MITO_SUBMIT_SCRIPTS_FROM} -le ${RAW_MITO_SUBMIT_SCRIPTS_TO} ]]
			then
				cd ../
				mkdir -p ${MITO_DIR}
				cd ${MITO_DIR}
				currentPhase=-1
				prefix=$(getPhaseFilePrefix)
				db=$(getCurrentDB)
        		currentStep=$((${RAW_MITO_SUBMIT_SCRIPTS_FROM}-1))
			elif [[ ${RAW_DASCOVER_SUBMIT_SCRIPTS_FROM} -gt 0 && ${RAW_DASCOVER_SUBMIT_SCRIPTS_FROM} -le ${RAW_DASCOVER_SUBMIT_SCRIPTS_TO} ]]
			then
				cd ../
				mkdir -p ${DASCOVER_DIR}
				cd ${DASCOVER_DIR}
				currentPhase=0
				prefix=$(getPhaseFilePrefix)
				db=$(getCurrentDB)
        		currentStep=$((${RAW_DASCOVER_SUBMIT_SCRIPTS_FROM}-1))				
			elif [[ ${RAW_REPMASK_SUBMIT_SCRIPTS_FROM} -gt 0 && ${RAW_REPMASK_SUBMIT_SCRIPTS_FROM} -le ${RAW_REPMASK_SUBMIT_SCRIPTS_TO} ]] ||
		   		[[ ${RAW_PATCH_SUBMIT_SCRIPTS_FROM} -gt 0 && ${RAW_PATCH_SUBMIT_SCRIPTS_FROM} -le ${RAW_PATCH_SUBMIT_SCRIPTS_TO} ]]
			then
				cd ../
				mkdir -p ${PATCHING_DIR}
				cd ${PATCHING_DIR}
				
				if [[ ${RAW_REPMASK_SUBMIT_SCRIPTS_FROM} -gt 0 && ${RAW_REPMASK_SUBMIT_SCRIPTS_FROM} -le ${RAW_REPMASK_SUBMIT_SCRIPTS_TO} ]]
				then
					currentPhase=1
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${RAW_REPMASK_SUBMIT_SCRIPTS_FROM}-1))
        		else
        			currentPhase=2
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${RAW_PATCH_SUBMIT_SCRIPTS_FROM}-1))
        		fi												
			elif [[ ${FIX_REPMASK_SUBMIT_SCRIPTS_FROM} -gt 0 && ${FIX_REPMASK_SUBMIT_SCRIPTS_FROM} -le ${FIX_REPMASK_SUBMIT_SCRIPTS_TO} ]] ||
		  	 	[[ ${FIX_SCRUB_SUBMIT_SCRIPTS_FROM} -gt 0 && ${FIX_SCRUB_SUBMIT_SCRIPTS_FROM} -le ${FIX_SCRUB_SUBMIT_SCRIPTS_TO} ]] ||
		   		[[ ${FIX_FILT_SUBMIT_SCRIPTS_FROM} -gt 0 && ${FIX_FILT_SUBMIT_SCRIPTS_FROM} -le ${FIX_FILT_SUBMIT_SCRIPTS_TO} ]] ||
		   		[[ ${FIX_TOUR_SUBMIT_SCRIPTS_FROM} -gt 0 && ${FIX_TOUR_SUBMIT_SCRIPTS_FROM} -le ${FIX_TOUR_SUBMIT_SCRIPTS_TO} ]] ||
		   		[[ ${FIX_CORR_SUBMIT_SCRIPTS_FROM} -gt 0 && ${FIX_CORR_SUBMIT_SCRIPTS_FROM} -le ${FIX_CORR_SUBMIT_SCRIPTS_TO} ]] ||
		   		[[ ${COR_CONTIG_SUBMIT_SCRIPTS_FROM} -gt 0 && ${COR_CONTIG_SUBMIT_SCRIPTS_FROM} -le ${COR_CONTIG_SUBMIT_SCRIPTS_TO} ]] ||
		   		[[ ${PB_ARROW_SUBMIT_SCRIPTS_FROM} -gt 0 && ${PB_ARROW_SUBMIT_SCRIPTS_FROM} -le ${PB_ARROW_SUBMIT_SCRIPTS_TO} ]] ||
		   		[[ ${CT_PURGEHAPLOTIGS_SUBMIT_SCRIPTS_FROM} -gt 0 && ${CT_PURGEHAPLOTIGS_SUBMIT_SCRIPTS_FROM} -le ${CT_PURGEHAPLOTIGS_SUBMIT_SCRIPTS_TO} ]] ||
		   		[[ ${CT_FREEBAYES_SUBMIT_SCRIPTS_FROM} -gt 0 && ${CT_FREEBAYES_SUBMIT_SCRIPTS_FROM} -le ${CT_FREEBAYES_SUBMIT_SCRIPTS_TO} ]] ||		   		
		   		[[ ${CT_PHASE_SUBMIT_SCRIPTS_FROM} -gt 0 && ${CT_PHASE_SUBMIT_SCRIPTS_FROM} -le ${CT_PHASE_SUBMIT_SCRIPTS_TO} ]] ||
		   		[[ ${SC_10X_SUBMIT_SCRIPTS_FROM} -gt 0 && ${SC_10X_SUBMIT_SCRIPTS_FROM} -le ${SC_10X_SUBMIT_SCRIPTS_TO} ]] ||
		   		[[ ${SC_BIONANO_SUBMIT_SCRIPTS_FROM} -gt 0 && ${SC_BIONANO_SUBMIT_SCRIPTS_FROM} -le ${SC_BIONANO_SUBMIT_SCRIPTS_TO} ]] ||
		   		[[ ${SC_HIC_SUBMIT_SCRIPTS_FROM} -gt 0 && ${SC_HIC_SUBMIT_SCRIPTS_FROM} -le ${SC_HIC_SUBMIT_SCRIPTS_TO} ]]		   		
			then
				cd ../	
			
		        if [[ -z "${FIX_REPMASK_USELAFIX_PATH}" ]]
				then 
					(>&2 echo "WARNING - Variable FIX_REPMASK_USELAFIX_PATH is not set.Try to use default path: patchedReads_dalign")
					FIX_REPMASK_USELAFIX_PATH="patchedReads_dalign"
				fi
				mkdir -p ${ASSMEBLY_DIR}/${FIX_REPMASK_USELAFIX_PATH}
				cd ${ASSMEBLY_DIR}/${FIX_REPMASK_USELAFIX_PATH}
				
				if [[ ${FIX_REPMASK_SUBMIT_SCRIPTS_FROM} -gt 0 && ${FIX_REPMASK_SUBMIT_SCRIPTS_FROM} -le ${FIX_REPMASK_SUBMIT_SCRIPTS_TO} ]]
				then
					currentPhase=3
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${FIX_REPMASK_SUBMIT_SCRIPTS_FROM}-1))
        		elif [[ ${FIX_SCRUB_SUBMIT_SCRIPTS_FROM} -gt 0 && ${FIX_SCRUB_SUBMIT_SCRIPTS_FROM} -le ${FIX_SCRUB_SUBMIT_SCRIPTS_TO} ]]
        		then
        			currentPhase=4
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${FIX_SCRUB_SUBMIT_SCRIPTS_FROM}-1))
        		elif [[ ${FIX_FILT_SUBMIT_SCRIPTS_FROM} -gt 0 && ${FIX_FILT_SUBMIT_SCRIPTS_FROM} -le ${FIX_FILT_SUBMIT_SCRIPTS_TO} ]]
        		then
        			currentPhase=5
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${FIX_FILT_SUBMIT_SCRIPTS_FROM}-1))
        		elif [[ ${FIX_TOUR_SUBMIT_SCRIPTS_FROM} -gt 0 && ${FIX_TOUR_SUBMIT_SCRIPTS_FROM} -le ${FIX_TOUR_SUBMIT_SCRIPTS_TO} ]]  
        		then
        			currentPhase=6
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${FIX_TOUR_SUBMIT_SCRIPTS_FROM}-1))
        		elif [[ ${FIX_CORR_SUBMIT_SCRIPTS_FROM} -gt 0 && ${FIX_CORR_SUBMIT_SCRIPTS_FROM} -le ${FIX_CORR_SUBMIT_SCRIPTS_TO} ]] 
        		then 
        			currentPhase=7
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${FIX_CORR_SUBMIT_SCRIPTS_FROM}-1))
        		elif ${COR_CONTIG_SUBMIT_SCRIPTS_FROM} -gt 0 && ${COR_CONTIG_SUBMIT_SCRIPTS_FROM} -le ${COR_CONTIG_SUBMIT_SCRIPTS_TO} ]]
        		then 
        			currentPhase=8
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${COR_CONTIG_SUBMIT_SCRIPTS_FROM}-1))
        		elif [[ ${PB_ARROW_SUBMIT_SCRIPTS_FROM} -gt 0 && ${PB_ARROW_SUBMIT_SCRIPTS_FROM} -le ${PB_ARROW_SUBMIT_SCRIPTS_TO} ]]
        		then
        			currentPhase=9
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${PB_ARROW_SUBMIT_SCRIPTS_FROM}-1))
        		elif [[ ${CT_PURGEHAPLOTIGS_SUBMIT_SCRIPTS_FROM} -gt 0 && ${CT_PURGEHAPLOTIGS_SUBMIT_SCRIPTS_FROM} -le ${CT_PURGEHAPLOTIGS_SUBMIT_SCRIPTS_TO} ]]
        		then 
        			currentPhase=10
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${CT_PURGEHAPLOTIGS_SUBMIT_SCRIPTS_FROM}-1))
        		elif [[ ${CT_FREEBAYES_SUBMIT_SCRIPTS_FROM} -gt 0 && ${CT_FREEBAYES_SUBMIT_SCRIPTS_FROM} -le ${CT_FREEBAYES_SUBMIT_SCRIPTS_TO} ]]
        		then 
        			currentPhase=11
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${CT_FREEBAYES_SUBMIT_SCRIPTS_FROM}-1))
        		elif [[ ${CT_PHASE_SUBMIT_SCRIPTS_FROM} -gt 0 && ${CT_PHASE_SUBMIT_SCRIPTS_FROM} -le ${CT_PHASE_SUBMIT_SCRIPTS_TO} ]]
        		then
        			currentPhase=12
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${CT_PHASE_SUBMIT_SCRIPTS_FROM}-1))
        		elif [[ ${SC_10X_SUBMIT_SCRIPTS_FROM} -gt 0 && ${SC_10X_SUBMIT_SCRIPTS_FROM} -le ${SC_10X_SUBMIT_SCRIPTS_TO} ]]
        		then
        			currentPhase=13
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${SC_10X_SUBMIT_SCRIPTS_FROM}-1))
        		elif [[ ${SC_BIONANO_SUBMIT_SCRIPTS_FROM} -gt 0 && ${SC_BIONANO_SUBMIT_SCRIPTS_FROM} -le ${SC_BIONANO_SUBMIT_SCRIPTS_TO} ]]
        		then
        			currentPhase=14
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${SC_BIONANO_SUBMIT_SCRIPTS_FROM}-1))
        		elif [[ ${SC_HIC_SUBMIT_SCRIPTS_FROM} -gt 0 && ${SC_HIC_SUBMIT_SCRIPTS_FROM} -le ${SC_HIC_SUBMIT_SCRIPTS_TO} ]]
        		then
        			currentPhase=15
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${SC_HIC_SUBMIT_SCRIPTS_FROM}-1))			
        		else
        			currentPhase=100 ## nothing to do, set phase to invalid value
				fi				
			fi 
		fi
	fi
fi

if [[ ${currentPhase} -eq -1 ]]
then
	if [[ $((${currentStep}+1)) -le ${RAW_MITO_SUBMIT_SCRIPTS_TO} ]]
    then
    	sbatch${appAccount} --job-name=${PROJECT_ID}_p${currentPhase}s$((${currentStep+1})) -o ${prefix}_step$((${currentStep}+1))_${db}_${slurmID}.out -e ${prefix}_step$((${currentStep}+1))_${db%.db}_${slurmID}.err -n1 -c1 -p ${SLURM_PARTITION} --time=01:00:00 --mem-per-cpu=6g --dependency=afterok:${RET##* } --wrap="bash ${SUBMIT_SCRIPTS_PATH}/createAndSubmitMarvelSlurmJobs.sh ${configFile} ${currentPhase} $((${currentStep}+1)) $slurmID"
        foundNext=1
	else
		### we have to create the coverage directory and change into that dir 
		if [[ $(echo "$(pwd)" | awk -F \/ '{print $NF}') -eq ${MITO_DIR} ]]
		then
			if [[ ${RAW_DASCOVER_SUBMIT_SCRIPTS_FROM} -gt 0 && ${RAW_DASCOVER_SUBMIT_SCRIPTS_FROM} -le ${RAW_DASCOVER_SUBMIT_SCRIPTS_TO} ]]
			then
				cd ../
				mkdir -p ${DASCOVER_DIR}
				cd ${DASCOVER_DIR}
				currentPhase=0
				prefix=$(getPhaseFilePrefix)
				db=$(getCurrentDB)
        		currentStep=$((${RAW_DASCOVER_SUBMIT_SCRIPTS_FROM}-1))				
			elif [[ ${RAW_REPMASK_SUBMIT_SCRIPTS_FROM} -gt 0 && ${RAW_REPMASK_SUBMIT_SCRIPTS_FROM} -le ${RAW_REPMASK_SUBMIT_SCRIPTS_TO} ]] ||
		   		[[ ${RAW_PATCH_SUBMIT_SCRIPTS_FROM} -gt 0 && ${RAW_PATCH_SUBMIT_SCRIPTS_FROM} -le ${RAW_PATCH_SUBMIT_SCRIPTS_TO} ]]
			then
				cd ../
				mkdir -p ${PATCHING_DIR}
				cd ${PATCHING_DIR}
				
				if [[ ${RAW_REPMASK_SUBMIT_SCRIPTS_FROM} -gt 0 && ${RAW_REPMASK_SUBMIT_SCRIPTS_FROM} -le ${RAW_REPMASK_SUBMIT_SCRIPTS_TO} ]]
				then
					currentPhase=1
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${RAW_REPMASK_SUBMIT_SCRIPTS_FROM}-1))
        		else
        			currentPhase=2
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${RAW_PATCH_SUBMIT_SCRIPTS_FROM}-1))
        		fi												
			elif [[ ${FIX_REPMASK_SUBMIT_SCRIPTS_FROM} -gt 0 && ${FIX_REPMASK_SUBMIT_SCRIPTS_FROM} -le ${FIX_REPMASK_SUBMIT_SCRIPTS_TO} ]] ||
		  	 	[[ ${FIX_SCRUB_SUBMIT_SCRIPTS_FROM} -gt 0 && ${FIX_SCRUB_SUBMIT_SCRIPTS_FROM} -le ${FIX_SCRUB_SUBMIT_SCRIPTS_TO} ]] ||
		   		[[ ${FIX_FILT_SUBMIT_SCRIPTS_FROM} -gt 0 && ${FIX_FILT_SUBMIT_SCRIPTS_FROM} -le ${FIX_FILT_SUBMIT_SCRIPTS_TO} ]] ||
		   		[[ ${FIX_TOUR_SUBMIT_SCRIPTS_FROM} -gt 0 && ${FIX_TOUR_SUBMIT_SCRIPTS_FROM} -le ${FIX_TOUR_SUBMIT_SCRIPTS_TO} ]] ||
		   		[[ ${FIX_CORR_SUBMIT_SCRIPTS_FROM} -gt 0 && ${FIX_CORR_SUBMIT_SCRIPTS_FROM} -le ${FIX_CORR_SUBMIT_SCRIPTS_TO} ]] ||
		   		[[ ${COR_CONTIG_SUBMIT_SCRIPTS_FROM} -gt 0 && ${COR_CONTIG_SUBMIT_SCRIPTS_FROM} -le ${COR_CONTIG_SUBMIT_SCRIPTS_TO} ]] ||
		   		[[ ${PB_ARROW_SUBMIT_SCRIPTS_FROM} -gt 0 && ${PB_ARROW_SUBMIT_SCRIPTS_FROM} -le ${PB_ARROW_SUBMIT_SCRIPTS_TO} ]] ||
		   		[[ ${CT_PURGEHAPLOTIGS_SUBMIT_SCRIPTS_FROM} -gt 0 && ${CT_PURGEHAPLOTIGS_SUBMIT_SCRIPTS_FROM} -le ${CT_PURGEHAPLOTIGS_SUBMIT_SCRIPTS_TO} ]] ||
		   		[[ ${CT_FREEBAYES_SUBMIT_SCRIPTS_FROM} -gt 0 && ${CT_FREEBAYES_SUBMIT_SCRIPTS_FROM} -le ${CT_FREEBAYES_SUBMIT_SCRIPTS_TO} ]] ||		   		
		   		[[ ${CT_PHASE_SUBMIT_SCRIPTS_FROM} -gt 0 && ${CT_PHASE_SUBMIT_SCRIPTS_FROM} -le ${CT_PHASE_SUBMIT_SCRIPTS_TO} ]] ||
		   		[[ ${SC_10X_SUBMIT_SCRIPTS_FROM} -gt 0 && ${SC_10X_SUBMIT_SCRIPTS_FROM} -le ${SC_10X_SUBMIT_SCRIPTS_TO} ]] ||
		   		[[ ${SC_BIONANO_SUBMIT_SCRIPTS_FROM} -gt 0 && ${SC_BIONANO_SUBMIT_SCRIPTS_FROM} -le ${SC_BIONANO_SUBMIT_SCRIPTS_TO} ]] ||
		   		[[ ${SC_HIC_SUBMIT_SCRIPTS_FROM} -gt 0 && ${SC_HIC_SUBMIT_SCRIPTS_FROM} -le ${SC_HIC_SUBMIT_SCRIPTS_TO} ]]		   		
			then
				cd ../	
			
		        if [[ -z "${FIX_REPMASK_USELAFIX_PATH}" ]]
				then 
					(>&2 echo "WARNING - Variable FIX_REPMASK_USELAFIX_PATH is not set.Try to use default path: patchedReads_dalign")
					FIX_REPMASK_USELAFIX_PATH="patchedReads_dalign"
				fi
				mkdir -p ${ASSMEBLY_DIR}/${FIX_REPMASK_USELAFIX_PATH}
				cd ${ASSMEBLY_DIR}/${FIX_REPMASK_USELAFIX_PATH}
				
				if [[ ${FIX_REPMASK_SUBMIT_SCRIPTS_FROM} -gt 0 && ${FIX_REPMASK_SUBMIT_SCRIPTS_FROM} -le ${FIX_REPMASK_SUBMIT_SCRIPTS_TO} ]]
				then
					currentPhase=3
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${FIX_REPMASK_SUBMIT_SCRIPTS_FROM}-1))
        		elif [[ ${FIX_SCRUB_SUBMIT_SCRIPTS_FROM} -gt 0 && ${FIX_SCRUB_SUBMIT_SCRIPTS_FROM} -le ${FIX_SCRUB_SUBMIT_SCRIPTS_TO} ]]
        		then
        			currentPhase=4
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${FIX_SCRUB_SUBMIT_SCRIPTS_FROM}-1))
        		elif [[ ${FIX_FILT_SUBMIT_SCRIPTS_FROM} -gt 0 && ${FIX_FILT_SUBMIT_SCRIPTS_FROM} -le ${FIX_FILT_SUBMIT_SCRIPTS_TO} ]]
        		then
        			currentPhase=5
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${FIX_FILT_SUBMIT_SCRIPTS_FROM}-1))
        		elif [[ ${FIX_TOUR_SUBMIT_SCRIPTS_FROM} -gt 0 && ${FIX_TOUR_SUBMIT_SCRIPTS_FROM} -le ${FIX_TOUR_SUBMIT_SCRIPTS_TO} ]]  
        		then
        			currentPhase=6
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${FIX_TOUR_SUBMIT_SCRIPTS_FROM}-1))
        		elif [[ ${FIX_CORR_SUBMIT_SCRIPTS_FROM} -gt 0 && ${FIX_CORR_SUBMIT_SCRIPTS_FROM} -le ${FIX_CORR_SUBMIT_SCRIPTS_TO} ]] 
        		then 
        			currentPhase=7
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${FIX_CORR_SUBMIT_SCRIPTS_FROM}-1))
        		elif ${COR_CONTIG_SUBMIT_SCRIPTS_FROM} -gt 0 && ${COR_CONTIG_SUBMIT_SCRIPTS_FROM} -le ${COR_CONTIG_SUBMIT_SCRIPTS_TO} ]]
        		then 
        			currentPhase=8
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${COR_CONTIG_SUBMIT_SCRIPTS_FROM}-1))
        		elif [[ ${PB_ARROW_SUBMIT_SCRIPTS_FROM} -gt 0 && ${PB_ARROW_SUBMIT_SCRIPTS_FROM} -le ${PB_ARROW_SUBMIT_SCRIPTS_TO} ]]
        		then
        			currentPhase=9
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${PB_ARROW_SUBMIT_SCRIPTS_FROM}-1))
        		elif [[ ${CT_PURGEHAPLOTIGS_SUBMIT_SCRIPTS_FROM} -gt 0 && ${CT_PURGEHAPLOTIGS_SUBMIT_SCRIPTS_FROM} -le ${CT_PURGEHAPLOTIGS_SUBMIT_SCRIPTS_TO} ]]
        		then 
        			currentPhase=10
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${CT_PURGEHAPLOTIGS_SUBMIT_SCRIPTS_FROM}-1))
        		elif [[ ${CT_FREEBAYES_SUBMIT_SCRIPTS_FROM} -gt 0 && ${CT_FREEBAYES_SUBMIT_SCRIPTS_FROM} -le ${CT_FREEBAYES_SUBMIT_SCRIPTS_TO} ]]
        		then 
        			currentPhase=11
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${CT_FREEBAYES_SUBMIT_SCRIPTS_FROM}-1))
        		elif [[ ${CT_PHASE_SUBMIT_SCRIPTS_FROM} -gt 0 && ${CT_PHASE_SUBMIT_SCRIPTS_FROM} -le ${CT_PHASE_SUBMIT_SCRIPTS_TO} ]]
        		then
        			currentPhase=12
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${CT_PHASE_SUBMIT_SCRIPTS_FROM}-1))
        		elif [[ ${SC_10X_SUBMIT_SCRIPTS_FROM} -gt 0 && ${SC_10X_SUBMIT_SCRIPTS_FROM} -le ${SC_10X_SUBMIT_SCRIPTS_TO} ]]
        		then
        			currentPhase=13
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${SC_10X_SUBMIT_SCRIPTS_FROM}-1))
        		elif [[ ${SC_BIONANO_SUBMIT_SCRIPTS_FROM} -gt 0 && ${SC_BIONANO_SUBMIT_SCRIPTS_FROM} -le ${SC_BIONANO_SUBMIT_SCRIPTS_TO} ]]
        		then
        			currentPhase=14
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${SC_BIONANO_SUBMIT_SCRIPTS_FROM}-1))
        		elif [[ ${SC_HIC_SUBMIT_SCRIPTS_FROM} -gt 0 && ${SC_HIC_SUBMIT_SCRIPTS_FROM} -le ${SC_HIC_SUBMIT_SCRIPTS_TO} ]]
        		then
        			currentPhase=15
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${SC_HIC_SUBMIT_SCRIPTS_FROM}-1))			
        		else
        			currentPhase=100 ## nothing to do, set phase to invalid value
				fi				
			fi 
		fi
	fi
fi

if [[ ${currentPhase} -eq 0 && ${foundNext} -eq 0 ]]
then 
    if [[ $((${currentStep}+1)) -le ${RAW_DASCOVER_SUBMIT_SCRIPTS_TO} ]]
    then 
        sbatch${appAccount} --job-name=${PROJECT_ID}_p${currentPhase}s$((${currentStep+1})) -o ${prefix}_step$((${currentStep}+1))_${RAW_DB%.db}_${slurmID}.out -e ${prefix}_step$((${currentStep}+1))_${RAW_DB%.db}_${slurmID}.err -n1 -c1 -p ${SLURM_PARTITION} --time=01:00:00 --mem-per-cpu=6g --dependency=afterok:${RET##* } --wrap="bash ${SUBMIT_SCRIPTS_PATH}/createAndSubmitMarvelSlurmJobs.sh ${configFile} ${currentPhase} $((${currentStep}+1)) $slurmID"
        foundNext=1
    else 
        currentPhase=1
        currentStep=$((${RAW_REPMASK_SUBMIT_SCRIPTS_FROM}-1))
		### we have to create the patching directory and change into that dir 
		if [[ $(echo "$(pwd)" | awk -F \/ '{print $NF}') -eq ${MITO_DIR} ]] || [[ $(echo "$(pwd)" | awk -F \/ '{print $NF}') -eq ${COVERAGE_DIR} ]]
		then
			if [[ ${RAW_REPMASK_SUBMIT_SCRIPTS_FROM} -gt 0 && ${RAW_REPMASK_SUBMIT_SCRIPTS_FROM} -le ${RAW_REPMASK_SUBMIT_SCRIPTS_TO} ]] ||
		   		[[ ${RAW_PATCH_SUBMIT_SCRIPTS_FROM} -gt 0 && ${RAW_PATCH_SUBMIT_SCRIPTS_FROM} -le ${RAW_PATCH_SUBMIT_SCRIPTS_TO} ]]
			then
				cd ../
				mkdir -p ${PATCHING_DIR}
				cd ${PATCHING_DIR}
				
				if [[ ${RAW_REPMASK_SUBMIT_SCRIPTS_FROM} -gt 0 && ${RAW_REPMASK_SUBMIT_SCRIPTS_FROM} -le ${RAW_REPMASK_SUBMIT_SCRIPTS_TO} ]]
				then
					currentPhase=1
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${RAW_REPMASK_SUBMIT_SCRIPTS_FROM}-1))
        		else
        			currentPhase=2
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${RAW_PATCH_SUBMIT_SCRIPTS_FROM}-1))
        		fi												
			elif [[ ${FIX_REPMASK_SUBMIT_SCRIPTS_FROM} -gt 0 && ${FIX_REPMASK_SUBMIT_SCRIPTS_FROM} -le ${FIX_REPMASK_SUBMIT_SCRIPTS_TO} ]] ||
		  	 	[[ ${FIX_SCRUB_SUBMIT_SCRIPTS_FROM} -gt 0 && ${FIX_SCRUB_SUBMIT_SCRIPTS_FROM} -le ${FIX_SCRUB_SUBMIT_SCRIPTS_TO} ]] ||
		   		[[ ${FIX_FILT_SUBMIT_SCRIPTS_FROM} -gt 0 && ${FIX_FILT_SUBMIT_SCRIPTS_FROM} -le ${FIX_FILT_SUBMIT_SCRIPTS_TO} ]] ||
		   		[[ ${FIX_TOUR_SUBMIT_SCRIPTS_FROM} -gt 0 && ${FIX_TOUR_SUBMIT_SCRIPTS_FROM} -le ${FIX_TOUR_SUBMIT_SCRIPTS_TO} ]] ||
		   		[[ ${FIX_CORR_SUBMIT_SCRIPTS_FROM} -gt 0 && ${FIX_CORR_SUBMIT_SCRIPTS_FROM} -le ${FIX_CORR_SUBMIT_SCRIPTS_TO} ]] ||
		   		[[ ${COR_CONTIG_SUBMIT_SCRIPTS_FROM} -gt 0 && ${COR_CONTIG_SUBMIT_SCRIPTS_FROM} -le ${COR_CONTIG_SUBMIT_SCRIPTS_TO} ]] ||
		   		[[ ${PB_ARROW_SUBMIT_SCRIPTS_FROM} -gt 0 && ${PB_ARROW_SUBMIT_SCRIPTS_FROM} -le ${PB_ARROW_SUBMIT_SCRIPTS_TO} ]] ||
		   		[[ ${CT_PURGEHAPLOTIGS_SUBMIT_SCRIPTS_FROM} -gt 0 && ${CT_PURGEHAPLOTIGS_SUBMIT_SCRIPTS_FROM} -le ${CT_PURGEHAPLOTIGS_SUBMIT_SCRIPTS_TO} ]] ||
		   		[[ ${CT_FREEBAYES_SUBMIT_SCRIPTS_FROM} -gt 0 && ${CT_FREEBAYES_SUBMIT_SCRIPTS_FROM} -le ${CT_FREEBAYES_SUBMIT_SCRIPTS_TO} ]] ||		   		
		   		[[ ${CT_PHASE_SUBMIT_SCRIPTS_FROM} -gt 0 && ${CT_PHASE_SUBMIT_SCRIPTS_FROM} -le ${CT_PHASE_SUBMIT_SCRIPTS_TO} ]] ||
		   		[[ ${SC_10X_SUBMIT_SCRIPTS_FROM} -gt 0 && ${SC_10X_SUBMIT_SCRIPTS_FROM} -le ${SC_10X_SUBMIT_SCRIPTS_TO} ]] ||
		   		[[ ${SC_BIONANO_SUBMIT_SCRIPTS_FROM} -gt 0 && ${SC_BIONANO_SUBMIT_SCRIPTS_FROM} -le ${SC_BIONANO_SUBMIT_SCRIPTS_TO} ]] ||
		   		[[ ${SC_HIC_SUBMIT_SCRIPTS_FROM} -gt 0 && ${SC_HIC_SUBMIT_SCRIPTS_FROM} -le ${SC_HIC_SUBMIT_SCRIPTS_TO} ]]		   		
			then
				cd ../	
			
		        if [[ -z "${FIX_REPMASK_USELAFIX_PATH}" ]]
				then 
					(>&2 echo "WARNING - Variable FIX_REPMASK_USELAFIX_PATH is not set.Try to use default path: patchedReads_dalign")
					FIX_REPMASK_USELAFIX_PATH="patchedReads_dalign"
				fi
				mkdir -p ${ASSMEBLY_DIR}/${FIX_REPMASK_USELAFIX_PATH}
				cd ${ASSMEBLY_DIR}/${FIX_REPMASK_USELAFIX_PATH}
				
				if [[ ${FIX_REPMASK_SUBMIT_SCRIPTS_FROM} -gt 0 && ${FIX_REPMASK_SUBMIT_SCRIPTS_FROM} -le ${FIX_REPMASK_SUBMIT_SCRIPTS_TO} ]]
				then
					currentPhase=3
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${FIX_REPMASK_SUBMIT_SCRIPTS_FROM}-1))
        		elif [[ ${FIX_SCRUB_SUBMIT_SCRIPTS_FROM} -gt 0 && ${FIX_SCRUB_SUBMIT_SCRIPTS_FROM} -le ${FIX_SCRUB_SUBMIT_SCRIPTS_TO} ]]
        		then
        			currentPhase=4
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${FIX_SCRUB_SUBMIT_SCRIPTS_FROM}-1))
        		elif [[ ${FIX_FILT_SUBMIT_SCRIPTS_FROM} -gt 0 && ${FIX_FILT_SUBMIT_SCRIPTS_FROM} -le ${FIX_FILT_SUBMIT_SCRIPTS_TO} ]]
        		then
        			currentPhase=5
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${FIX_FILT_SUBMIT_SCRIPTS_FROM}-1))
        		elif [[ ${FIX_TOUR_SUBMIT_SCRIPTS_FROM} -gt 0 && ${FIX_TOUR_SUBMIT_SCRIPTS_FROM} -le ${FIX_TOUR_SUBMIT_SCRIPTS_TO} ]]  
        		then
        			currentPhase=6
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${FIX_TOUR_SUBMIT_SCRIPTS_FROM}-1))
        		elif [[ ${FIX_CORR_SUBMIT_SCRIPTS_FROM} -gt 0 && ${FIX_CORR_SUBMIT_SCRIPTS_FROM} -le ${FIX_CORR_SUBMIT_SCRIPTS_TO} ]] 
        		then 
        			currentPhase=7
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${FIX_CORR_SUBMIT_SCRIPTS_FROM}-1))
        		elif ${COR_CONTIG_SUBMIT_SCRIPTS_FROM} -gt 0 && ${COR_CONTIG_SUBMIT_SCRIPTS_FROM} -le ${COR_CONTIG_SUBMIT_SCRIPTS_TO} ]]
        		then 
        			currentPhase=8
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${COR_CONTIG_SUBMIT_SCRIPTS_FROM}-1))
        		elif [[ ${PB_ARROW_SUBMIT_SCRIPTS_FROM} -gt 0 && ${PB_ARROW_SUBMIT_SCRIPTS_FROM} -le ${PB_ARROW_SUBMIT_SCRIPTS_TO} ]]
        		then
        			currentPhase=9
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${PB_ARROW_SUBMIT_SCRIPTS_FROM}-1))
        		elif [[ ${CT_PURGEHAPLOTIGS_SUBMIT_SCRIPTS_FROM} -gt 0 && ${CT_PURGEHAPLOTIGS_SUBMIT_SCRIPTS_FROM} -le ${CT_PURGEHAPLOTIGS_SUBMIT_SCRIPTS_TO} ]]
        		then 
        			currentPhase=10
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${CT_PURGEHAPLOTIGS_SUBMIT_SCRIPTS_FROM}-1))
        		elif [[ ${CT_FREEBAYES_SUBMIT_SCRIPTS_FROM} -gt 0 && ${CT_FREEBAYES_SUBMIT_SCRIPTS_FROM} -le ${CT_FREEBAYES_SUBMIT_SCRIPTS_TO} ]]
        		then 
        			currentPhase=11
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${CT_FREEBAYES_SUBMIT_SCRIPTS_FROM}-1))
        		elif [[ ${CT_PHASE_SUBMIT_SCRIPTS_FROM} -gt 0 && ${CT_PHASE_SUBMIT_SCRIPTS_FROM} -le ${CT_PHASE_SUBMIT_SCRIPTS_TO} ]]
        		then
        			currentPhase=12
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${CT_PHASE_SUBMIT_SCRIPTS_FROM}-1))
				elif [[ ${SC_10X_SUBMIT_SCRIPTS_FROM} -gt 0 && ${SC_10X_SUBMIT_SCRIPTS_FROM} -le ${SC_10X_SUBMIT_SCRIPTS_TO} ]]
        		then
        			currentPhase=13
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${SC_10X_SUBMIT_SCRIPTS_FROM}-1))  
        		elif [[ ${SC_BIONANO_SUBMIT_SCRIPTS_FROM} -gt 0 && ${SC_BIONANO_SUBMIT_SCRIPTS_FROM} -le ${SC_BIONANO_SUBMIT_SCRIPTS_TO} ]]
        		then
        			currentPhase=14
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${SC_BIONANO_SUBMIT_SCRIPTS_FROM}-1))
        		elif [[ ${SC_HIC_SUBMIT_SCRIPTS_FROM} -gt 0 && ${SC_HIC_SUBMIT_SCRIPTS_FROM} -le ${SC_HIC_SUBMIT_SCRIPTS_TO} ]]
        		then
        			currentPhase=15
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${SC_HIC_SUBMIT_SCRIPTS_FROM}-1))        				      			
        		else
        			currentPhase=100 ## nothing to do, set phase to invalid value
				fi				
			fi
		fi
    fi
fi

if [[ ${currentPhase} -eq 1  && ${foundNext} -eq 0 ]]
then 
    if [[ $((${currentStep}+1)) -gt 0 && $((${currentStep}+1)) -le ${RAW_REPMASK_SUBMIT_SCRIPTS_TO} ]]
    then 
        sbatch${appAccount} --job-name=${PROJECT_ID}_p${currentPhase}s$((${currentStep+1})) -o ${prefix}_step$((${currentStep}+1))_${RAW_DB%.db}_${slurmID}.out -e ${prefix}_step$((${currentStep}+1))_${RAW_DB%.db}_${slurmID}.err -n1 -c1 -p ${SLURM_PARTITION} --time=01:00:00 --mem-per-cpu=6g --dependency=afterok:${RET##* } --wrap="bash ${SUBMIT_SCRIPTS_PATH}/createAndSubmitMarvelSlurmJobs.sh ${configFile} ${currentPhase} $((${currentStep}+1)) $slurmID"
        foundNext=1
    else 
        currentPhase=2
        currentStep=$((${RAW_PATCH_SUBMIT_SCRIPTS_FROM}-1))
    fi
fi  

if [[ ${currentPhase} -eq 2 && ${foundNext} -eq 0 ]]
then 
    if [[ $((${currentStep}+1)) -gt 0 && $((${currentStep}+1)) -le ${RAW_PATCH_SUBMIT_SCRIPTS_TO} ]]
    then 
        sbatch${appAccount} --job-name=${PROJECT_ID}_p${currentPhase}s$((${currentStep+1})) -o ${prefix}_step$((${currentStep}+1))_${RAW_DB%.db}_${slurmID}.out -e ${prefix}_step$((${currentStep}+1))_${RAW_DB%.db}_${slurmID}.err -n1 -c1 -p ${SLURM_PARTITION} --time=01:00:00 --mem-per-cpu=6g --dependency=afterok:${RET##* } --wrap="bash ${SUBMIT_SCRIPTS_PATH}/createAndSubmitMarvelSlurmJobs.sh ${configFile} ${currentPhase} $((${currentStep}+1)) $slurmID"
        foundNext=1
    else 
        currentPhase=3
        currentStep=$((${FIX_REPMASK_SUBMIT_SCRIPTS_FROM}-1))
        ### we have to create the assembly directory and change into that dir 
		if [[ $(echo "$(pwd)" | awk -F \/ '{print $NF}') -eq ${MITO_DIR} ]] || [[ $(echo "$(pwd)" | awk -F \/ '{print $NF}') -eq ${COVERAGE_DIR} ]] || [[ $(echo "$(pwd)" | awk -F \/ '{print $NF}') -eq ${PATCHING_DIR} ]]
		then
			if [[ ${FIX_REPMASK_SUBMIT_SCRIPTS_FROM} -gt 0 && ${FIX_REPMASK_SUBMIT_SCRIPTS_FROM} -le ${FIX_REPMASK_SUBMIT_SCRIPTS_TO} ]] ||
		  	 	[[ ${FIX_SCRUB_SUBMIT_SCRIPTS_FROM} -gt 0 && ${FIX_SCRUB_SUBMIT_SCRIPTS_FROM} -le ${FIX_SCRUB_SUBMIT_SCRIPTS_TO} ]] ||
		   		[[ ${FIX_FILT_SUBMIT_SCRIPTS_FROM} -gt 0 && ${FIX_FILT_SUBMIT_SCRIPTS_FROM} -le ${FIX_FILT_SUBMIT_SCRIPTS_TO} ]] ||
		   		[[ ${FIX_TOUR_SUBMIT_SCRIPTS_FROM} -gt 0 && ${FIX_TOUR_SUBMIT_SCRIPTS_FROM} -le ${FIX_TOUR_SUBMIT_SCRIPTS_TO} ]] ||
		   		[[ ${FIX_CORR_SUBMIT_SCRIPTS_FROM} -gt 0 && ${FIX_CORR_SUBMIT_SCRIPTS_FROM} -le ${FIX_CORR_SUBMIT_SCRIPTS_TO} ]] ||
		   		[[ ${COR_CONTIG_SUBMIT_SCRIPTS_FROM} -gt 0 && ${COR_CONTIG_SUBMIT_SCRIPTS_FROM} -le ${COR_CONTIG_SUBMIT_SCRIPTS_TO} ]] ||
		   		[[ ${PB_ARROW_SUBMIT_SCRIPTS_FROM} -gt 0 && ${PB_ARROW_SUBMIT_SCRIPTS_FROM} -le ${PB_ARROW_SUBMIT_SCRIPTS_TO} ]] ||
		   		[[ ${CT_PURGEHAPLOTIGS_SUBMIT_SCRIPTS_FROM} -gt 0 && ${CT_PURGEHAPLOTIGS_SUBMIT_SCRIPTS_FROM} -le ${CT_PURGEHAPLOTIGS_SUBMIT_SCRIPTS_TO} ]] ||
		   		[[ ${CT_FREEBAYES_SUBMIT_SCRIPTS_FROM} -gt 0 && ${CT_FREEBAYES_SUBMIT_SCRIPTS_FROM} -le ${CT_FREEBAYES_SUBMIT_SCRIPTS_TO} ]] ||		   		
		   		[[ ${CT_PHASE_SUBMIT_SCRIPTS_FROM} -gt 0 && ${CT_PHASE_SUBMIT_SCRIPTS_FROM} -le ${CT_PHASE_SUBMIT_SCRIPTS_TO} ]] ||
		   		[[ ${SC_10X_SUBMIT_SCRIPTS_FROM} -gt 0 && ${SC_10X_SUBMIT_SCRIPTS_FROM} -le ${SC_10X_SUBMIT_SCRIPTS_TO} ]] ||
		   		[[ ${SC_BIONANO_SUBMIT_SCRIPTS_FROM} -gt 0 && ${SC_BIONANO_SUBMIT_SCRIPTS_FROM} -le ${SC_BIONANO_SUBMIT_SCRIPTS_TO} ]] ||
		   		[[ ${SC_HIC_SUBMIT_SCRIPTS_FROM} -gt 0 && ${SC_HIC_SUBMIT_SCRIPTS_FROM} -le ${SC_HIC_SUBMIT_SCRIPTS_TO} ]]		   		
			then
				cd ../	
			
		        if [[ -z "${FIX_REPMASK_USELAFIX_PATH}" ]]
				then 
					(>&2 echo "WARNING - Variable FIX_REPMASK_USELAFIX_PATH is not set.Try to use default path: patchedReads_dalign")
					FIX_REPMASK_USELAFIX_PATH="patchedReads_dalign"
				fi
				mkdir -p ${ASSMEBLY_DIR}/${FIX_REPMASK_USELAFIX_PATH}
				cd ${ASSMEBLY_DIR}/${FIX_REPMASK_USELAFIX_PATH}
				
				if [[ ${FIX_REPMASK_SUBMIT_SCRIPTS_FROM} -gt 0 && ${FIX_REPMASK_SUBMIT_SCRIPTS_FROM} -le ${FIX_REPMASK_SUBMIT_SCRIPTS_TO} ]]
				then
					currentPhase=3
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${FIX_REPMASK_SUBMIT_SCRIPTS_FROM}-1))
        		elif [[ ${FIX_SCRUB_SUBMIT_SCRIPTS_FROM} -gt 0 && ${FIX_SCRUB_SUBMIT_SCRIPTS_FROM} -le ${FIX_SCRUB_SUBMIT_SCRIPTS_TO} ]]
        		then
        			currentPhase=4
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${FIX_SCRUB_SUBMIT_SCRIPTS_FROM}-1))
        		elif [[ ${FIX_FILT_SUBMIT_SCRIPTS_FROM} -gt 0 && ${FIX_FILT_SUBMIT_SCRIPTS_FROM} -le ${FIX_FILT_SUBMIT_SCRIPTS_TO} ]]
        		then
        			currentPhase=5
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${FIX_FILT_SUBMIT_SCRIPTS_FROM}-1))
        		elif [[ ${FIX_TOUR_SUBMIT_SCRIPTS_FROM} -gt 0 && ${FIX_TOUR_SUBMIT_SCRIPTS_FROM} -le ${FIX_TOUR_SUBMIT_SCRIPTS_TO} ]]  
        		then
        			currentPhase=6
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${FIX_TOUR_SUBMIT_SCRIPTS_FROM}-1))
        		elif [[ ${FIX_CORR_SUBMIT_SCRIPTS_FROM} -gt 0 && ${FIX_CORR_SUBMIT_SCRIPTS_FROM} -le ${FIX_CORR_SUBMIT_SCRIPTS_TO} ]] 
        		then 
        			currentPhase=7
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${FIX_CORR_SUBMIT_SCRIPTS_FROM}-1))
        		elif ${COR_CONTIG_SUBMIT_SCRIPTS_FROM} -gt 0 && ${COR_CONTIG_SUBMIT_SCRIPTS_FROM} -le ${COR_CONTIG_SUBMIT_SCRIPTS_TO} ]]
        		then 
        			currentPhase=8
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${COR_CONTIG_SUBMIT_SCRIPTS_FROM}-1))
        		elif [[ ${PB_ARROW_SUBMIT_SCRIPTS_FROM} -gt 0 && ${PB_ARROW_SUBMIT_SCRIPTS_FROM} -le ${PB_ARROW_SUBMIT_SCRIPTS_TO} ]]
        		then
        			currentPhase=9
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${PB_ARROW_SUBMIT_SCRIPTS_FROM}-1))
        		elif [[ ${CT_PURGEHAPLOTIGS_SUBMIT_SCRIPTS_FROM} -gt 0 && ${CT_PURGEHAPLOTIGS_SUBMIT_SCRIPTS_FROM} -le ${CT_PURGEHAPLOTIGS_SUBMIT_SCRIPTS_TO} ]]
        		then 
        			currentPhase=10
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${CT_PURGEHAPLOTIGS_SUBMIT_SCRIPTS_FROM}-1))
        		elif [[ ${CT_FREEBAYES_SUBMIT_SCRIPTS_FROM} -gt 0 && ${CT_FREEBAYES_SUBMIT_SCRIPTS_FROM} -le ${CT_FREEBAYES_SUBMIT_SCRIPTS_TO} ]]
        		then 
        			currentPhase=11
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${CT_FREEBAYES_SUBMIT_SCRIPTS_FROM}-1))
        		elif [[ ${CT_PHASE_SUBMIT_SCRIPTS_FROM} -gt 0 && ${CT_PHASE_SUBMIT_SCRIPTS_FROM} -le ${CT_PHASE_SUBMIT_SCRIPTS_TO} ]]
        		then
        			currentPhase=12
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${CT_PHASE_SUBMIT_SCRIPTS_FROM}-1))
        		elif [[ ${SC_10X_SUBMIT_SCRIPTS_FROM} -gt 0 && ${SC_10X_SUBMIT_SCRIPTS_FROM} -le ${SC_10X_SUBMIT_SCRIPTS_TO} ]]
        		then
        			currentPhase=13
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${SC_10X_SUBMIT_SCRIPTS_FROM}-1))
        		elif [[ ${SC_BIONANO_SUBMIT_SCRIPTS_FROM} -gt 0 && ${SC_BIONANO_SUBMIT_SCRIPTS_FROM} -le ${SC_BIONANO_SUBMIT_SCRIPTS_TO} ]]
        		then
        			currentPhase=14
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${SC_BIONANO_SUBMIT_SCRIPTS_FROM}-1))
        		elif [[ ${SC_HIC_SUBMIT_SCRIPTS_FROM} -gt 0 && ${SC_HIC_SUBMIT_SCRIPTS_FROM} -le ${SC_HIC_SUBMIT_SCRIPTS_TO} ]]
        		then
        			currentPhase=15
					prefix=$(getPhaseFilePrefix)
					db=$(getCurrentDB)
        			currentStep=$((${SC_HIC_SUBMIT_SCRIPTS_FROM}-1))		
        		else
        			currentPhase=100 ## nothing to do, set phase to invalid value
				fi				
			fi
		fi
    fi 
fi  

if [[ ${currentPhase} -eq 3 && ${foundNext} -eq 0 ]]
then 
    if [[ $((${currentStep}+1)) -gt 0 && $((${currentStep}+1)) -le ${FIX_REPMASK_SUBMIT_SCRIPTS_TO} ]]
    then 
        sbatch${appAccount} --job-name=${PROJECT_ID}_p${currentPhase}s$((${currentStep+1})) -o ${prefix}_step$((${currentStep}+1))_${RAW_DB%.db}_${slurmID}.out -e ${prefix}_step$((${currentStep}+1))_${RAW_DB%.db}_${slurmID}.err -n1 -c1 -p ${SLURM_PARTITION} --time=01:00:00 --mem-per-cpu=6g --dependency=afterok:${RET##* } --wrap="bash ${SUBMIT_SCRIPTS_PATH}/createAndSubmitMarvelSlurmJobs.sh ${configFile} ${currentPhase} $((${currentStep}+1)) $slurmID"
        foundNext=1
    else 
        currentPhase=4
        currentStep=$((${FIX_SCRUB_SUBMIT_SCRIPTS_FROM}-1))
    fi
fi

if [[ ${currentPhase} -eq 4 && ${foundNext} -eq 0 ]]
then 
    if [[ $((${currentStep}+1)) -gt 0 && $((${currentStep}+1)) -le ${FIX_SCRUB_SUBMIT_SCRIPTS_TO} ]]
    then 
        sbatch${appAccount} --job-name=${PROJECT_ID}_p${currentPhase}s$((${currentStep+1})) -o ${prefix}_step$((${currentStep}+1))_${FIX_DB%.db}_${slurmID}.out -e ${prefix}_step$((${currentStep}+1))_${FIX_DB%.db}_${slurmID}.err -n1 -c1 -p ${SLURM_PARTITION} --time=01:00:00 --mem-per-cpu=6g --dependency=afterok:${RET##* } --wrap="bash ${SUBMIT_SCRIPTS_PATH}/createAndSubmitMarvelSlurmJobs.sh ${configFile} ${currentPhase} $((${currentStep}+1)) $slurmID"
        foundNext=1
    else 
        currentPhase=5
        currentStep=$((${FIX_FILT_SUBMIT_SCRIPTS_FROM}-1))
    fi 
fi     

if [[ ${currentPhase} -eq 5 && ${foundNext} -eq 0 ]]
then 
    if [[ ${FIX_FILT_TYPE} -eq 0 && ${currentStep} -eq 2 && -n ${FIX_FILT_LAFILTER_RMSYMROUNDS} && ${FIX_FILT_LAFILTER_RMSYMROUNDS} -gt 0 && ! -f filt_02_LAfilter_block_${FIX_DB%.db}.${slurmID}.plan ]]
    then                 
        sbatch${appAccount} --job-name=${PROJECT_ID}_p${currentPhase}s${currentStep} -o ${prefix}_step$((${currentStep}+1))_${FIX_DB%.db}_${slurmID}.out -e ${prefix}_step$((${currentStep}+1))_${FIX_DB%.db}_${slurmID}.err -n1 -c1 -p ${SLURM_PARTITION} --time=01:00:00 --mem-per-cpu=6g --dependency=afterok:${RET##* } --wrap="bash ${SUBMIT_SCRIPTS_PATH}/createAndSubmitMarvelSlurmJobs.sh ${configFile} ${currentPhase} ${currentStep} $slurmID"
        foundNext=1
    elif [[ $((${currentStep}+1)) -gt 0 && $((${currentStep}+1)) -le ${FIX_FILT_SUBMIT_SCRIPTS_TO} ]]
    then 
        sbatch${appAccount} --job-name=${PROJECT_ID}_p${currentPhase}s$((${currentStep+1})) -o ${prefix}_step$((${currentStep}+1))_${FIX_DB%.db}_${slurmID}.out -e ${prefix}_step$((${currentStep}+1))_${FIX_DB%.db}_${slurmID}.err -n1 -c1 -p ${SLURM_PARTITION} --time=01:00:00 --mem-per-cpu=6g --dependency=afterok:${RET##* } --wrap="bash ${SUBMIT_SCRIPTS_PATH}/createAndSubmitMarvelSlurmJobs.sh ${configFile} ${currentPhase} $((${currentStep}+1)) $slurmID"
        foundNext=1
    else 
        currentPhase=6
        currentStep=$((${FIX_TOUR_SUBMIT_SCRIPTS_FROM}-1))
    fi 
fi

if [[ ${currentPhase} -eq 6 && ${foundNext} -eq 0 ]]
then 
    if [[ $((${currentStep}+1)) -gt 0 && $((${currentStep}+1)) -le ${FIX_TOUR_SUBMIT_SCRIPTS_TO} ]]
    then 
        sbatch${appAccount} --job-name=${PROJECT_ID}_p${currentPhase}s$((${currentStep+1})) -o ${prefix}_step$((${currentStep}+1))_${FIX_DB%.db}_${slurmID}.out -e ${prefix}_step$((${currentStep}+1))_${FIX_DB%.db}_${slurmID}.err -n1 -c1 -p ${SLURM_PARTITION} --time=01:00:00 --mem-per-cpu=6g --dependency=afterok:${RET##* } --wrap="bash ${SUBMIT_SCRIPTS_PATH}/createAndSubmitMarvelSlurmJobs.sh ${configFile} ${currentPhase} $((${currentStep}+1)) $slurmID"
        foundNext=1
    else 
        currentPhase=7
        currentStep=$((${FIX_CORR_SUBMIT_SCRIPTS_FROM}-1))
    fi
fi 

if [[ ${currentPhase} -eq 7 && ${foundNext} -eq 0 ]]
then 
    if [[ $((${currentStep}+1)) -gt 0 && $((${currentStep}+1)) -le ${FIX_CORR_SUBMIT_SCRIPTS_TO} ]]
    then 
        sbatch${appAccount} --job-name=${PROJECT_ID}_p${currentPhase}s$((${currentStep+1})) -o ${prefix}_step$((${currentStep}+1))_${FIX_DB%.db}_${slurmID}.out -e ${prefix}_step$((${currentStep}+1))_${FIX_DB%.db}_${slurmID}.err -n1 -c1 -p ${SLURM_PARTITION} --time=01:00:00 --mem-per-cpu=6g --dependency=afterok:${RET##* } --wrap="bash ${SUBMIT_SCRIPTS_PATH}/createAndSubmitMarvelSlurmJobs.sh ${configFile} ${currentPhase} $((${currentStep}+1)) $slurmID"
        foundNext=1
    else 
        currentPhase=8
        currentStep=$((${COR_CONTIG_SUBMIT_SCRIPTS_FROM}-1))    
    fi
fi

if [[ ${currentPhase} -eq 8 && ${foundNext} -eq 0 ]]
then 
    if [[ $((${currentStep}+1)) -gt 0 && $((${currentStep}+1)) -le ${COR_CONTIG_SUBMIT_SCRIPTS_TO} ]]
    then 
        sbatch${appAccount} --job-name=${PROJECT_ID}_p${currentPhase}s$((${currentStep+1})) -o ${prefix}_step$((${currentStep}+1))_${CONT_DB%.db}_${slurmID}.out -e ${prefix}_step$((${currentStep}+1))_${CONT_DB%.db}_${slurmID}.err -n1 -c1 -p ${SLURM_PARTITION} --time=01:00:00 --mem-per-cpu=6g --dependency=afterok:${RET##* } --wrap="bash ${SUBMIT_SCRIPTS_PATH}/createAndSubmitMarvelSlurmJobs.sh ${configFile} ${currentPhase} $((${currentStep}+1)) $slurmID"
        foundNext=1
	else 
        currentPhase=9
        currentStep=$((${PB_ARROW_SUBMIT_SCRIPTS_FROM}-1))        
    fi
fi

if [[ ${currentPhase} -eq 9 && ${foundNext} -eq 0 ]]
then 
    if [[ $((${currentStep}+1)) -gt 0 && $((${currentStep}+1)) -le ${PB_ARROW_SUBMIT_SCRIPTS_TO} ]]
    then 
        sbatch${appAccount} --job-name=${PROJECT_ID}_p${currentPhase}s$((${currentStep+1})) -o ${prefix}_step$((${currentStep}+1))_${RAW_DB%.db}_${slurmID}.out -e ${prefix}_step$((${currentStep}+1))_${RAW_DB%.db}_${slurmID}.err -n1 -c1 -p ${SLURM_PARTITION} --time=01:00:00 --mem-per-cpu=6g --dependency=afterok:${RET##* } --wrap="bash ${SUBMIT_SCRIPTS_PATH}/createAndSubmitMarvelSlurmJobs.sh ${configFile} ${currentPhase} $((${currentStep}+1)) $slurmID"
        foundNext=1
	else 
        currentPhase=10
        currentStep=$((${CT_PURGEHAPLOTIGS_SUBMIT_SCRIPTS_FROM}-1))                
    fi
fi

if [[ ${currentPhase} -eq 10 && ${foundNext} -eq 0 ]]
then 
    if [[ $((${currentStep}+1)) -gt 0 && $((${currentStep}+1)) -le ${CT_PURGEHAPLOTIGS_SUBMIT_SCRIPTS_TO} ]]
    then 
        sbatch${appAccount} --job-name=${PROJECT_ID}_p${currentPhase}s$((${currentStep+1})) -o ${prefix}_step$((${currentStep}+1))_${FIX_DB%.db}_${slurmID}.out -e ${prefix}_step$((${currentStep}+1))_${FIX_DB%.db}_${slurmID}.err -n1 -c1 -p ${SLURM_PARTITION} --time=01:00:00 --mem-per-cpu=6g --dependency=afterok:${RET##* } --wrap="bash ${SUBMIT_SCRIPTS_PATH}/createAndSubmitMarvelSlurmJobs.sh ${configFile} ${currentPhase} $((${currentStep}+1)) $slurmID"
        foundNext=1
	else 
        currentPhase=11
        currentStep=$((${CT_FREEBAYES_SUBMIT_SCRIPTS_FROM}-1))                        
    fi
fi

if [[ ${currentPhase} -eq 11 && ${foundNext} -eq 0 ]]
then 
    if [[ $((${currentStep}+1)) -gt 0 && $((${currentStep}+1)) -le ${CT_FREEBAYES_SUBMIT_SCRIPTS_TO} ]]
    then 
        sbatch${appAccount} --job-name=${PROJECT_ID}_p${currentPhase}s$((${currentStep+1})) -o ${prefix}_step$((${currentStep}+1))_${FIX_DB%.db}_${slurmID}.out -e ${prefix}_step$((${currentStep}+1))_${FIX_DB%.db}_${slurmID}.err -n1 -c1 -p ${SLURM_PARTITION} --time=01:00:00 --mem-per-cpu=6g --dependency=afterok:${RET##* } --wrap="bash ${SUBMIT_SCRIPTS_PATH}/createAndSubmitMarvelSlurmJobs.sh ${configFile} ${currentPhase} $((${currentStep}+1)) $slurmID"
        foundNext=1
	else 
        currentPhase=12
        currentStep=$((${CT_PHASE_SUBMIT_SCRIPTS_FROM}-1))        
    fi
fi

if [[ ${currentPhase} -eq 12 && ${foundNext} -eq 0 ]]
then 
    if [[ $((${currentStep}+1)) -gt 0 && $((${currentStep}+1)) -le ${CT_PHASE_SUBMIT_SCRIPTS_TO} ]]
    then 
        sbatch${appAccount} --job-name=${PROJECT_ID}_p${currentPhase}s$((${currentStep+1})) -o ${prefix}_step$((${currentStep}+1))_${FIX_DB%.db}_${slurmID}.out -e ${prefix}_step$((${currentStep}+1))_${FIX_DB%.db}_${slurmID}.err -n1 -c1 -p ${SLURM_PARTITION} --time=01:00:00 --mem-per-cpu=6g --dependency=afterok:${RET##* } --wrap="bash ${SUBMIT_SCRIPTS_PATH}/createAndSubmitMarvelSlurmJobs.sh ${configFile} ${currentPhase} $((${currentStep}+1)) $slurmID"
        foundNext=1
    else 
        currentPhase=13
        currentStep=$((${SC_10X_SUBMIT_SCRIPTS_FROM}-1))    
    fi
fi                        

if [[ ${currentPhase} -eq 13 && ${foundNext} -eq 0 ]]
then 
    if [[ $((${currentStep}+1)) -gt 0 && $((${currentStep}+1)) -le ${SC_10X_SUBMIT_SCRIPTS_TO} ]]
    then 
        sbatch${appAccount} --job-name=${PROJECT_ID}_p${currentPhase}s$((${currentStep+1})) -o ${prefix}_step$((${currentStep}+1))_${FIX_DB%.db}_${slurmID}.out -e ${prefix}_step$((${currentStep}+1))_${FIX_DB%.db}_${slurmID}.err -n1 -c1 -p ${SLURM_PARTITION} --time=01:00:00 --mem-per-cpu=6g --dependency=afterok:${RET##* } --wrap="bash ${SUBMIT_SCRIPTS_PATH}/createAndSubmitMarvelSlurmJobs.sh ${configFile} ${currentPhase} $((${currentStep}+1)) $slurmID"
        foundNext=1
    else 
        currentPhase=14
        currentStep=$((${SC_BIONANO_SUBMIT_SCRIPTS_FROM}-1))    
    fi
fi

if [[ ${currentPhase} -eq 14 && ${foundNext} -eq 0 ]]
then 
    if [[ $((${currentStep}+1)) -gt 0 && $((${currentStep}+1)) -le ${SC_BIONANO_SUBMIT_SCRIPTS_TO} ]]
    then 
        sbatch${appAccount} --job-name=${PROJECT_ID}_p${currentPhase}s$((${currentStep+1})) -o ${prefix}_step$((${currentStep}+1))_${FIX_DB%.db}_${slurmID}.out -e ${prefix}_step$((${currentStep}+1))_${FIX_DB%.db}_${slurmID}.err -n1 -c1 -p ${SLURM_PARTITION} --time=01:00:00 --mem-per-cpu=6g --dependency=afterok:${RET##* } --wrap="bash ${SUBMIT_SCRIPTS_PATH}/createAndSubmitMarvelSlurmJobs.sh ${configFile} ${currentPhase} $((${currentStep}+1)) $slurmID"
        foundNext=1
    else 
        currentPhase=15
        currentStep=$((${SC_HIC_SUBMIT_SCRIPTS_FROM}-1))    
    fi
fi   

if [[ ${currentPhase} -eq 15 && ${foundNext} -eq 0 ]]
then 
    if [[ $((${currentStep}+1)) -gt 0 && $((${currentStep}+1)) -le ${SC_HIC_SUBMIT_SCRIPTS_TO} ]]
    then 
        sbatch${appAccount} --job-name=${PROJECT_ID}_p${currentPhase}s$((${currentStep+1})) -o ${prefix}_step$((${currentStep}+1))_${FIX_DB%.db}_${slurmID}.out -e ${prefix}_step$((${currentStep}+1))_${FIX_DB%.db}_${slurmID}.err -n1 -c1 -p ${SLURM_PARTITION} --time=01:00:00 --mem-per-cpu=6g --dependency=afterok:${RET##* } --wrap="bash ${SUBMIT_SCRIPTS_PATH}/createAndSubmitMarvelSlurmJobs.sh ${configFile} ${currentPhase} $((${currentStep}+1)) $slurmID"
        foundNext=1
    fi
fi                     
                        

if [[ ${foundNext} -eq 0 ]]
then
	# submit a dummy job that waits until the last real jobs sucessfully finished
sbatch${appAccount} --job-name=${PROJECT_ID}_final -o ${prefix}_final_step_${currentPhase}_${currentStep}_${slurmID}.out -e ${prefix}_final_step_${currentPhase}_${currentStep}_${slurmID}.err -n1 -c1 -p ${SLURM_PARTITION} --time=00:15:00 --mem=1g --dependency=afterok:${RET##* } --wrap="sleep 5 && echo \"finished - all selected jobs created and submitted. Last Step: ${prefix} ${currentPhase} ${currentStep} $slurmID ${configFile}\""     
fi