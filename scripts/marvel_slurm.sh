#!/bin/bash 

configFile=$1
Id=$2
if [[ -z "$Id" ]]
then
  Id=1
fi

if [[ ! -f ${configFile} ]]
then 
    (>&2 echo "cannot access config file ${configFile}")
    exit 1
fi

source ${configFile}

## do some general sanity checks
	
if [[ -z "${PROJECT_ID}" ]]
then 
    (>&2 echo "ERROR - You have to specify a project id. Set variable PROJECT_ID")
    exit 1
fi

if [[ -d "${PATCHING_DIR}" ]]
then 
    (>&2 echo "ERROR - You have to set PATCHING_DIR.")
    exit 1
fi

if [[ -z "${ASSMEBLY_DIR}" ]]
then 
    (>&2 echo "ERROR - You have to set ASSMEBLY_DIR")
    exit 1
fi

if [[ -z "${DB_PATH}" ]]
then 
    (>&2 echo "ERROR - You have to set DB_PATH. Location of the initial databases MARVEL and DAZZLER.")
    exit 1
fi

## find entry point to create first plan and submit that stuff 

if [[ ${RAW_REPMASK_SUBMIT_SCRIPTS_FROM} -gt 0 ]] 
then 
    currentPhase=0
    currentStep=${RAW_REPMASK_SUBMIT_SCRIPTS_FROM}    
elif [[ ${RAW_PATCH_SUBMIT_SCRIPTS_FROM} -gt 0 ]] 
then 
    currentPhase=1
    currentStep=${RAW_PATCH_SUBMIT_SCRIPTS_FROM}    
elif [[ ${FIX_SCRUB_SUBMIT_SCRIPTS_FROM} -gt 0 ]] 
then 
    currentPhase=2
    currentStep=${FIX_SCRUB_SUBMIT_SCRIPTS_FROM}        
elif [[ ${FIX_FILT_SUBMIT_SCRIPTS_FROM} -gt 0 ]] 
then 
    currentPhase=3
    currentStep=${FIX_FILT_SUBMIT_SCRIPTS_FROM}        
elif [[ ${FIX_TOUR_SUBMIT_SCRIPTS_FROM} -gt 0 ]] 
then 
    currentPhase=4
    currentStep=${FIX_TOUR_SUBMIT_SCRIPTS_FROM}        
elif [[ ${FIX_CORR_SUBMIT_SCRIPTS_FROM} -gt 0 ]] 
then 
    currentPhase=5
    currentStep=${FIX_CORR_SUBMIT_SCRIPTS_FROM}        
elif [[ ${COR_CONTIG_SUBMIT_SCRIPTS_FROM} -gt 0 ]] 
then 
    currentPhase=6
    currentStep=${COR_CONTIG_SUBMIT_SCRIPTS_FROM}        
elif [[ ${PB_ARROW_SUBMIT_SCRIPTS_FROM} -gt 0 ]] 
then 
    currentPhase=7
    currentStep=${PB_ARROW_SUBMIT_SCRIPTS_FROM}        
else 
    echo "nothing to do"
    exit 0
fi

${SUBMIT_SCRIPTS_PATH}/createAndSubmitMarvelSlurmJobs.sh ${configFile} ${currentPhase} ${currentStep} ${Id}
