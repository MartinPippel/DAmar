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
## phases 0-DAScover, 1-repmask, 2-patching, 3-repmask, 4-scrubbing, 5-filtering, 6-touring, 7-correction, 8-contigAnalysis, 9-arrow 

if [[ ${currentPhase} -eq -2 ]]
then	 
	${SUBMIT_SCRIPTS_PATH}/createQCandStatsPlans.sh ${configFile} ${currentStep} ${slurmID}
    if [ $? -ne 0 ]
    then 
        (>&2 echo "createQCandStatsPlans.sh failed some how. Stop here.")
        exit 1      
    fi
elif [[ ${currentPhase} -eq -1 ]]
then	 
	if [[ ! -f ${RAW_DB%db}.db ]]
	then 
		if [[ ! -f ${DB_PATH}/${RAW_DB%db}.db ]]
		then 
			(>&2 echo "Cannot find initial databases ${RAW_DB%db}.db in directory ${DB_PATH}")
	        exit 1	
		fi		
		cp ${DB_PATH}/${RAW_DB%db}.db ${DB_PATH}/.${RAW_DB%db}.idx ${DB_PATH}/.${RAW_DB%db}.bps .
		if [[ -f ${DB_PATH}/.${RAW_DB%db}.pacbio.anno && -f ${DB_PATH}/.${RAW_DB%db}.pacbio.data ]]
		then 
			cp 	${DB_PATH}/.${RAW_DB%db}.pacbio.anno ${DB_PATH}/.${RAW_DB%db}.pacbio.data .
		fi
		if [[ -f ${DB_PATH}/.${RAW_DB%db}.seqID.anno && -f ${DB_PATH}/.${RAW_DB%db}.seqID.data ]]
		then 
			cp 	${DB_PATH}/.${RAW_DB%db}.seqID.anno ${DB_PATH}/.${RAW_DB%db}.seqID.data .
		fi
	fi		
    ${SUBMIT_SCRIPTS_PATH}/createMitoAssemblyPlans.sh ${configFile} ${currentStep} ${slurmID}
    if [ $? -ne 0 ]
    then 
        (>&2 echo "createMitoAssemblyPlans.sh failed some how. Stop here.")
        exit 1      
    fi
elif [[ ${currentPhase} -eq 0 ]]
then	 
	if [[ ! -f ${RAW_DAZZ_DB%db}.db ]]
	then 
		if [[ ! -f ${DB_PATH}/${RAW_DB%db}.db || ! -f ${DB_PATH}/${RAW_DAZZ_DB%db}.db ]]
		then 
			(>&2 echo "Cannot find initial databases ${RAW_DB%db}.db and ${RAW_DAZZ_DB%db}.db in directory ${DB_PATH}")
	        exit 1	
		fi		
		ln -s -r ${DB_PATH}/${RAW_DAZZ_DB%db}.db ${DB_PATH}/.${RAW_DAZZ_DB%db}.idx ${DB_PATH}/.${RAW_DAZZ_DB%db}.bps .
	fi		
    ${SUBMIT_SCRIPTS_PATH}/createDAScoverPlans.sh ${configFile} ${currentStep} ${slurmID}
    if [ $? -ne 0 ]
    then 
        (>&2 echo "createDAScoverPlans.sh failed some how. Stop here.")
        exit 1      
    fi
elif [[ ${currentPhase} -eq 1 ]]
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
	if [[ ${currentStep} -eq 1 ]]
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
		if [[ ! -d "${FIX_REPMASK_USELAFIX_PATH}" ]]	
		then 
			if [[ ! -d "../${PATCHING_DIR}/${FIX_REPMASK_USELAFIX_PATH}" ]]
			then
				(>&2 echo "Cannot find patched reads in directory ../${PATCHING_DIR}/${FIX_REPMASK_USELAFIX_PATH}")
				(>&2 echo "cwd $(pwd)")
		        exit 1
			fi	
			ln -s -r ../${PATCHING_DIR}/${FIX_REPMASK_USELAFIX_PATH} ${FIX_REPMASK_USELAFIX_PATH}
		fi
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
elif [[ ${currentPhase} -eq 10 ]]
then 
    ${SUBMIT_SCRIPTS_PATH}/createPurgeHaplotigPlans.sh ${configFile} ${currentStep} ${slurmID}
    if [ $? -ne 0 ]
    then 
        (>&2 echo "${SUBMIT_SCRIPTS_PATH}/createPurgeHaplotigPlans.sh failed some how. Stop here.")
        exit 1      
    fi      
elif [[ ${currentPhase} -eq 11 ]]
then 
    ${SUBMIT_SCRIPTS_PATH}/createFreeBayesPolishingPlans.sh ${configFile} ${currentStep} ${slurmID}
    if [ $? -ne 0 ]
    then 
        (>&2 echo "${SUBMIT_SCRIPTS_PATH}/createFreeBayesPolishingPlans.sh failed some how. Stop here.")
        exit 1      
    fi
elif [[ ${currentPhase} -eq 12 ]]
then 
    ${SUBMIT_SCRIPTS_PATH}/createPhasePlans.sh ${configFile} ${currentStep} ${slurmID}
    if [ $? -ne 0 ]
    then 
        (>&2 echo "${SUBMIT_SCRIPTS_PATH}/createPhasePlans.sh failed some how. Stop here.")
        exit 1      
    fi
elif [[ ${currentPhase} -eq 13 ]]
then 
    ${SUBMIT_SCRIPTS_PATH}/create10XPlans.sh ${configFile} ${currentStep} ${slurmID}
    if [ $? -ne 0 ]
    then 
        (>&2 echo "${SUBMIT_SCRIPTS_PATH}/createScaff10XPlans.sh failed some how. Stop here.")
        exit 1      
    fi      
elif [[ ${currentPhase} -eq 14 ]]
then 
    ${SUBMIT_SCRIPTS_PATH}/createBionanoPlans.sh ${configFile} ${currentStep} ${slurmID}
    if [ $? -ne 0 ]
    then 
        (>&2 echo "${SUBMIT_SCRIPTS_PATH}/createBionanoPlans.sh failed some how. Stop here.")
        exit 1      
    fi      
elif [[ ${currentPhase} -eq 15 ]]
then 
    ${SUBMIT_SCRIPTS_PATH}/createHiCPlans.sh ${configFile} ${currentStep} ${slurmID}
    if [ $? -ne 0 ]
    then 
        (>&2 echo "${SUBMIT_SCRIPTS_PATH}/createHiCPlans.sh failed some how. Stop here.")
        exit 1      
    fi          
else
    (>&2 echo "unknown assembly phase: ${currentPhase}")
    exit 1
fi 