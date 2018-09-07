#!/bin/bash 

configFile=$1
currentPhase=$2
currentStep=$3
slurmID=$4

cwd=$(pwd)
echo "createCommandPlan.sh ${configFile} ${currentPhase} ${currentStep} ${slurmID}"
echo "cwd ${cwd}" 

if [[ ! -f ${configFile} ]]
then 
    (>&2 echo "cannot access config file ${configFile}")
    exit 1
fi

source ${configFile}

## todo sanity checks 
## phases 0-repmask, 1-patching, 2-scrubbing, 3-filtering, 4-touring, 5-correction 


if [[ ${currentPhase} -eq 1 ]]
then	 
	if [[ ! -f ${RAW_DB%db}.db ]]
	then 
		if [[ ! -f ${DB_PATH}/${RAW_DB%db}.db || ! -f ${DB_PATH}/${RAW_DAZZ_DB%db}.db ]]
		then 
			(>&2 echo "Cannot find initial databases ${RAW_DB%db}.db and ${RAW_DAZZ_DB%db}.db in directory ${DB_PATH}")
	        exit 1	
		fi		
		ln -s -r ${DB_PATH}/${RAW_DB%db}.db ${DB_PATH}/.${RAW_DB%db}.idx ${DB_PATH}/.${RAW_DB%db}.bps .
		ln -s -r ${DB_PATH}/${RAW_DAZZ_DB%db}.db ${DB_PATH}/.${RAW_DAZZ_DB%db}.idx ${DB_PATH}/.${RAW_DAZZ_DB%db}.bps .
	fi		
    ${SUBMIT_SCRIPTS_PATH}/createRepmaskPlans.sh ${configFile} ${currentStep} ${slurmID}
    if [ $? -ne 0 ]
    then 
        (>&2 echo "createRepmaskPlans.sh failed some how. Stop here.")
        exit 1      
    fi
elif [[ ${currentPhase} -eq 2 ]]
then 
	${SUBMIT_SCRIPTS_PATH}/createReadPatchingPlans.sh ${configFile} ${currentStep} ${slurmID}
    if [ $? -ne 0 ]
    then 
        (>&2 echo "createReadPatchingPlans.sh failed some how. Stop here.")
        exit 1      
    fi
elif [[ ${currentPhase} -eq 3 ]]    
then 
	if [[ -z ${RAW_FIX_LAFIX_PATH} ]]
	then 
		(>&2 echo "Variable RAW_FIX_LAFIX_PATH must be set. Its used to create the assembly subfolder ASSMEBLY_DIR/RAW_FIX_LAFIX_PATH")
	    exit 1
	fi
	if [[ ! -f ${RAW_DB%db}.db ]]
	then 
		if [[ ! -f ${DB_PATH}/${RAW_DB%db}.db || ! -f ${DB_PATH}/${RAW_DAZZ_DB%db}.db ]]
		then 
			(>&2 echo "Cannot find initial databases ${RAW_DB%db}.db and ${RAW_DAZZ_DB%db}.db in directory ${DB_PATH}")
	        exit 1	
		fi		
		ln -s -r ${DB_PATH}/${RAW_DB%db}.db ${DB_PATH}/.${RAW_DB%db}.idx ${DB_PATH}/.${RAW_DB%db}.bps .
		ln -s -r ${DB_PATH}/${RAW_DAZZ_DB%db}.db ${DB_PATH}/.${RAW_DAZZ_DB%db}.idx ${DB_PATH}/.${RAW_DAZZ_DB%db}.bps .
	fi
    ${SUBMIT_SCRIPTS_PATH}/createRepmaskPlans2.sh ${configFile} ${currentStep} ${slurmID}
    if [ $? -ne 0 ]
    then 
        (>&2 echo "createScrubbingPlans2.sh failed some how. Stop here.")
        exit 1      
    fi   
elif [[ ${currentPhase} -eq 4 ]]    
then 
    ${SUBMIT_SCRIPTS_PATH}/createScrubbingPlans.sh ${configFile} ${currentStep} ${slurmID}
    if [ $? -ne 0 ]
    then 
        (>&2 echo "createScrubbingPlans.sh failed some how. Stop here.")
        exit 1      
    fi
elif [[ ${currentPhase} -eq 5 ]]
then 
    ${SUBMIT_SCRIPTS_PATH}/createFilteringPlans.sh ${configFile} ${currentStep} ${slurmID}
    if [ $? -ne 0 ]
    then 
        (>&2 echo "createFilteringPlans.sh failed some how. Stop here.")
        exit 1      
    fi
elif [[ ${currentPhase} -eq 6 ]]
then 
    ${SUBMIT_SCRIPTS_PATH}/createTouringPlans.sh ${configFile} ${currentStep} ${slurmID}
    if [ $? -ne 0 ]
    then 
        (>&2 echo "createTouringPlans.sh failed some how. Stop here.")
        exit 1      
    fi
elif [[ ${currentPhase} -eq 7 ]]
then 
	${SUBMIT_SCRIPTS_PATH}/createCorrectionPlans.sh ${configFile} ${currentStep} ${slurmID}
    if [ $? -ne 0 ]
    then 
        (>&2 echo "createCorrectionPlans.sh failed some how. Stop here.")
        exit 1      
    fi
elif [[ ${currentPhase} -eq 8 ]]
then 
	${SUBMIT_SCRIPTS_PATH}/createContigAnalyzePlans.sh ${configFile} ${currentStep} ${slurmID}
    if [ $? -ne 0 ]
    then 
        (>&2 echo "${SUBMIT_SCRIPTS_PATH}/createContigAnalyzePlans.sh failed some how. Stop here.")
        exit 1      
    fi
elif [[ ${currentPhase} -eq 9 ]]
then 
    ${SUBMIT_SCRIPTS_PATH}/createPacBioArrowPlans.sh ${configFile} ${currentStep} ${slurmID}
    if [ $? -ne 0 ]
    then 
        (>&2 echo "${SUBMIT_SCRIPTS_PATH}/createPacBioArrowPlans.sh failed some how. Stop here.")
        exit 1      
    fi 
else
    (>&2 echo "unknown assembly phase: ${currentPhase}")
    exit 1
fi 