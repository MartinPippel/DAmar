#!/bin/bash 

config=$1
curJobID=$2
lastJobID=$3
slurmID=$4

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

source ${config}

  c=${curJobID}
  cjobid=$(prependLeadingNull $c)
  if ! ls scrub_${cjobid}_*_*_${FIX_DB%.db}.${slurmID}.plan 1> /dev/null 2>&1;
  then
   echo "missing file: scrub_${cjobid}_*_*_${FIX_DB%.db}.${slurmID}.plan" 
   exit 1
  fi
  # get job name 
  jname=$(ls scrub_${cjobid}_*_*_${FIX_DB%.db}.${slurmID}.plan | awk -F \_ '{print $3}')
  jtype=$(ls scrub_${cjobid}_*_*_${FIX_DB%.db}.${slurmID}.plan | awk -F \_ '{print $4}')
  if [[ ${jtype} != "block" && ${jtype} != "single" ]]
  then
    echo "unknown slurm type: ${jtype}. valid types: block, single"
    exit 1
  fi

  REP_MEM=${MEM_DEFAULT}
  TMP="MEM_${jname}"
  if [[ -n ${!TMP} ]]
  then 
    REP_MEM=${!TMP}
  fi
  REP_TIME=${TIME_DEFAULT}
  TMP="TIME_${jname}"
  if [[ -n ${!TMP} ]]
  then
    REP_TIME=${!TMP}
  fi
  REP_CORES=${THREADS_DEFAULT}
  TMP="THREADS_${jname}"
  if [[ -n ${!TMP} ]]
  then
    REP_CORES=${!TMP}
  fi
  REP_NTASKS_PER_NODE=""
  TMP="TASKS_${jname}"
  if [[ -n ${!TMP} ]]
  then
    REP_NTASKS_PER_NODE="#SBATCH --ntasks-per-node=${!TMP}"
  fi
  STEPSIZE=1
  TMP="STEPSIZE_${jname}"
  if [[ -n ${!TMP} ]]
  then
    if [[ ${!TMP} -gt 1 ]]
    then
        STEPSIZE="${!TMP}"
    fi
  fi

  REP_JOBS=$(wc -l scrub_${cjobid}_${jname}_${jtype}_${FIX_DB%.db}.${slurmID}.plan | awk '{print $1}')
  mkdir -p log_scrub_${jname}_${FIX_DB%.db}  
  first=1
  if [[ ${REP_JOBS} -gt 10000 ]]
  then
    from=1
    to=10000
    d=1
    while [ $from -lt ${REP_JOBS} ]
    do
      sed -n ${from},${to}p scrub_${cjobid}_${jname}_${jtype}_${FIX_DB%.db}.plan > scrub_${cjobid}_${jname}_${jtype}_${FIX_DB%.db}.${d}.${slurmID}.plan
      jobs=$((${to}-${from}+1))
      if [[ ${STEPSIZE} -gt 1 ]]
      then
        jobs="${jobs}:${STEPSIZE}"
      fi
    file=scrub_${cjobid}_${jname}_${jtype}_${FIX_DB%.db}.${d}.${slurmID}.slurm
sed -e "s:#REP_NTASKS_PER_NODE:${REP_NTASKS_PER_NODE}:" \
-e "s:REP_NAME:${FIX_DB%.db}:g" \
-e "s/REP_TIME/${REP_TIME}/" \
-e "s/REP_MEMORY/${REP_MEM}/g" \
-e "s/REP_JOBS/${jobs}/" \
-e "s:REP_JNAME:scrub_${cjobid}_${FIX_DB%.db}_${d}:" \
-e "s:REP_CORES:${REP_CORES}:" \
-e "s:REP_SCRIPT:scrub_${cjobid}_${jname}_${jtype}_${FIX_DB%.db}.${d}.${slurmID}.plan:" \
-e "s:REP_PARTITION:${FIX_SLURM_PARTITION}:" \
-e "s:REP_LOG:log_scrub_${jname}_${FIX_DB%.db}:" ${SUBMIT_SCRIPTS_PATH}/${jtype}.slurm > ${file}

      if [[ ${first} -eq 1 ]]
      then
        first=0
        RET=$(sbatch ${file}) && isNumber ${RET##* } && echo "submit ${file} ${RET##* }"        
      else
        RET=$(sbatch --dependency=afterok:${RET##* } ${file}) && isNumber ${RET##* } echo "submit ${file} ${RET##* }"
      fi
      d=$(($d+1))
      from=$((${to}+1))
      to=$((${to}+10000))
      if [[ $to -gt ${REP_JOBS} ]]
      then
        to=${REP_JOBS}
      fi
      sleep 5
    done
  else
      if [[ ${STEPSIZE} -gt 1 ]]
      then
        REP_JOBS="${REP_JOBS}:${STEPSIZE}"
      fi
file=scrub_${cjobid}_${jname}_${jtype}_${FIX_DB%.db}.${slurmID}.slurm
sed -e "s:#REP_NTASKS_PER_NODE:${REP_NTASKS_PER_NODE}:" \
-e "s:REP_NAME:${FIX_DB%.db}:g" \
-e "s/REP_TIME/${REP_TIME}/" \
-e "s/REP_MEMORY/${REP_MEM}/g" \
-e "s/REP_JOBS/${REP_JOBS}/" \
-e "s:REP_JNAME:scrub_${cjobid}_${FIX_DB%.db}:" \
-e "s:REP_CORES:${REP_CORES}:" \
-e "s:REP_SCRIPT:scrub_${cjobid}_${jname}_${jtype}_${FIX_DB%.db}.${slurmID}.plan:" \
-e "s:REP_PARTITION:${FIX_SLURM_PARTITION}:" \
-e "s:REP_LOG:log_scrub_${jname}_${FIX_DB%.db}:" ${SUBMIT_SCRIPTS_PATH}/${jtype}.slurm > $file
RET=$(sbatch ${file}) && isNumber ${RET##* } && echo "submit ${file} ${RET##* }"
fi


c=$((c+1))
if [[ $c -le $lastJobID ]]
then
    sbatch --job-name=scrub_step${c}_${FIX_DB%.db} -o scrub_step${c}_${FIX_DB%.db}.out -e scrub_step${c}_${FIX_DB%.db}.err -n1 -c1 -p batch --time=01:00:00 --mem=12g --dependency=afterok:${RET} --wrap="sh ${SUBMIT_SCRIPTS_PATH}/createSubmitScripts.sh ${config} $c $lastJobID $slurmID" 
fi 