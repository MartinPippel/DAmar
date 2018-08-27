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
else 
    echo "nothing to do"
    exit 0
fi

${SUBMIT_SCRIPTS_PATH}/createAndSubmitMarvelSlurmJobs.sh ${configFile} ${currentPhase} ${currentStep} ${Id}
