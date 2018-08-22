#!/bin/bash 

configFile=$1
currentPhase=$2
currentStep=$3
slurmID=$4

echo "createCommandPlan.sh ${configFile} ${currentPhase} ${currentStep} ${slurmID}"
echo $(pwd) 

if [[ ! -f ${configFile} ]]
then 
    (>&2 echo "cannot access config file ${configFile}")
    exit 1
fi

source ${configFile}

## todo sanity checks 
## phases 0-repmask, 1-patching, 2-scrubbing, 3-filtering, 4-touring, 5-correction 

if [[ ${currentPhase} -eq 0 ]]
then 
    ${SUBMIT_SCRIPTS_PATH}/createRepmaskPlans.sh ${configFile} ${currentStep} ${slurmID}
    if [ $? -ne 0 ]
    then 
        (>&2 echo "createRepmaskPlans.sh failed some how. Stop here.")
        exit 1      
    fi
elif [[ ${currentPhase} -eq 1 ]]
then 
    ${SUBMIT_SCRIPTS_PATH}/createReadPatchingPlans.sh ${configFile} ${currentStep} ${slurmID}
    if [ $? -ne 0 ]
    then 
        (>&2 echo "createReadPatchingPlans.sh failed some how. Stop here.")
        exit 1      
    fi
elif [[ ${currentPhase} -eq 2 ]]    
then 
    ${SUBMIT_SCRIPTS_PATH}/createScrubbingPlans.sh ${configFile} ${currentStep} ${slurmID}
    if [ $? -ne 0 ]
    then 
        (>&2 echo "createScrubbingPlans.sh failed some how. Stop here.")
        exit 1      
    fi
elif [[ ${currentPhase} -eq 3 ]]
then 
    ${SUBMIT_SCRIPTS_PATH}/createFilteringPlans.sh ${configFile} ${currentStep} ${slurmID}
    if [ $? -ne 0 ]
    then 
        (>&2 echo "createFilteringPlans.sh failed some how. Stop here.")
        exit 1      
    fi
elif [[ ${currentPhase} -eq 4 ]]
then 
    ${SUBMIT_SCRIPTS_PATH}/createTouringPlans.sh ${configFile} ${currentStep} ${slurmID}
    if [ $? -ne 0 ]
    then 
        (>&2 echo "createTouringPlans.sh failed some how. Stop here.")
        exit 1      
    fi
elif [[ ${currentPhase} -eq 5 ]]
then 
    ${SUBMIT_SCRIPTS_PATH}/createCorrectionPlans.sh ${configFile} ${currentStep} ${slurmID}
    if [ $? -ne 0 ]
    then 
        (>&2 echo "createCorrectionPlans.sh failed some how. Stop here.")
        exit 1      
    fi
else
    echo "unknown assembly phase: ${currentPhase}"
    exit 1
fi 