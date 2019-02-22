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

if [[ -z "${COVERAGE_DIR}" ]]
then 
    (>&2 echo "ERROR - You have to set COVERAGE_DIR.")
    exit 1
fi

if [[ -z "${PATCHING_DIR}" ]]
then 
    (>&2 echo "ERROR - You have to set PATCHING_DIR.")
    exit 1
fi

if [[ -z "${ASSMEBLY_DIR}" ]]
then 
    (>&2 echo "ERROR - You have to set ASSMEBLY_DIR")
    exit 1
fi

if [[ "${ASSMEBLY_DIR}" == "${PATCHING_DIR}" ]]
then 
    (>&2 echo "ERROR - PATCHING_DIR must be different from ASSMEBLY_DIR")
    exit 1
fi

if [[ -z "${DB_PATH}" ]]
then 
    (>&2 echo "ERROR - You have to set DB_PATH. Location of the initial databases MARVEL and DAZZLER.")
    exit 1
fi

if [[ -z "${QC_DATA_DIR}" ]]
then 
    (>&2 echo "ERROR - You have to set QC_DATA_DIR")
    exit 1
fi

# file must be present 
function realpath()
{
	echo "$(cd "$(dirname "$1")"; pwd)/$(basename "$1")"
}

## find entry point to create first plan and submit that stuff 

if [[ ${RAW_QC_SUBMIT_SCRIPTS_FROM} -gt 0 ]] 
then 
    currentPhase=-2
    currentStep=${RAW_QC_SUBMIT_SCRIPTS_FROM}     
elif [[ ${RAW_MITO_SUBMIT_SCRIPTS_FROM} -gt 0 ]] 
then 
    currentPhase=-1
    currentStep=${RAW_MITO_SUBMIT_SCRIPTS_FROM}    
elif [[ ${RAW_DASCOVER_SUBMIT_SCRIPTS_FROM} -gt 0 ]] 
then 
    currentPhase=0
    currentStep=${RAW_DASCOVER_SUBMIT_SCRIPTS_FROM}    
elif [[ ${RAW_REPMASK_SUBMIT_SCRIPTS_FROM} -gt 0 ]] 
then 
    currentPhase=1
    currentStep=${RAW_REPMASK_SUBMIT_SCRIPTS_FROM}    
elif [[ ${RAW_PATCH_SUBMIT_SCRIPTS_FROM} -gt 0 ]] 
then 
    currentPhase=2
    currentStep=${RAW_PATCH_SUBMIT_SCRIPTS_FROM} 
elif [[ ${FIX_REPMASK_SUBMIT_SCRIPTS_FROM} -gt 0 ]] 
then 
    currentPhase=3
    currentStep=${FIX_REPMASK_SUBMIT_SCRIPTS_FROM}    
    
elif [[ ${FIX_SCRUB_SUBMIT_SCRIPTS_FROM} -gt 0 ]] 
then 
    currentPhase=4
    currentStep=${FIX_SCRUB_SUBMIT_SCRIPTS_FROM}        
elif [[ ${FIX_FILT_SUBMIT_SCRIPTS_FROM} -gt 0 ]] 
then 
    currentPhase=5
    currentStep=${FIX_FILT_SUBMIT_SCRIPTS_FROM}        
elif [[ ${FIX_TOUR_SUBMIT_SCRIPTS_FROM} -gt 0 ]] 
then 
    currentPhase=6
    currentStep=${FIX_TOUR_SUBMIT_SCRIPTS_FROM}        
elif [[ ${FIX_CORR_SUBMIT_SCRIPTS_FROM} -gt 0 ]] 
then 
    currentPhase=7
    currentStep=${FIX_CORR_SUBMIT_SCRIPTS_FROM}        
elif [[ ${COR_CONTIG_SUBMIT_SCRIPTS_FROM} -gt 0 ]] 
then 
    currentPhase=8
    currentStep=${COR_CONTIG_SUBMIT_SCRIPTS_FROM}        
elif [[ ${PB_ARROW_SUBMIT_SCRIPTS_FROM} -gt 0 ]] 
then 
    currentPhase=9
    currentStep=${PB_ARROW_SUBMIT_SCRIPTS_FROM}
elif [[ ${CT_PURGEHAPLOTIGS_SUBMIT_SCRIPTS_FROM} -gt 0 ]] 
then 
    currentPhase=10
    currentStep=${CT_PURGEHAPLOTIGS_SUBMIT_SCRIPTS_FROM}
elif [[ ${CT_FREEBAYES_SUBMIT_SCRIPTS_FROM} -gt 0 ]] 
then 
    currentPhase=11
    currentStep=${CT_FREEBAYES_SUBMIT_SCRIPTS_FROM}                                           
elif [[ ${CT_WHATSHAP_SUBMIT_SCRIPTS_FROM} -gt 0 ]] 
then 
    currentPhase=12
    currentStep=${CT_WHATSHAP_SUBMIT_SCRIPTS_FROM}
elif [[ ${SC_10X_SUBMIT_SCRIPTS_FROM} -gt 0 ]] 
then 
    currentPhase=13
    currentStep=${SC_10X_SUBMIT_SCRIPTS_FROM}
elif [[ ${SC_BIONANO_SUBMIT_SCRIPTS_FROM} -gt 0 ]] 
then 
    currentPhase=14
    currentStep=${SC_BIONANO_SUBMIT_SCRIPTS_FROM}
elif [[ ${SC_HIC_SUBMIT_SCRIPTS_FROM} -gt 0 ]] 
then 
    currentPhase=15
    currentStep=${SC_HIC_SUBMIT_SCRIPTS_FROM}        
else 
    echo "nothing to do"
    exit 0
fi

realPathConfigFile=$(realpath "${configFile}")

cwd=$(pwd)

if [[ ${currentPhase} -eq -2 ]]
then 
	mkdir -p ${QC_DATA_DIR}
	cd ${QC_DATA_DIR}
	${SUBMIT_SCRIPTS_PATH}/createAndSubmitMarvelSlurmJobs.sh ${realPathConfigFile} ${currentPhase} ${currentStep} ${Id}
	cd ${cwd}
elif [[ ${currentPhase} -eq -1 ]]
then 
	mkdir -p ${MITO_DIR}
	cd ${MITO_DIR}
	${SUBMIT_SCRIPTS_PATH}/createAndSubmitMarvelSlurmJobs.sh ${realPathConfigFile} ${currentPhase} ${currentStep} ${Id}
	cd ${cwd}
elif [[ ${currentPhase} -eq 0 ]]
then 
	mkdir -p ${COVERAGE_DIR}
	cd ${COVERAGE_DIR}
	${SUBMIT_SCRIPTS_PATH}/createAndSubmitMarvelSlurmJobs.sh ${realPathConfigFile} ${currentPhase} ${currentStep} ${Id}
	cd ${cwd}
elif [[ ${currentPhase} -lt 3 ]]
then 
	mkdir -p ${PATCHING_DIR}
	cd ${PATCHING_DIR}
	${SUBMIT_SCRIPTS_PATH}/createAndSubmitMarvelSlurmJobs.sh ${realPathConfigFile} ${currentPhase} ${currentStep} ${Id}
	cd ${cwd}
elif [[ ${currentPhase} -lt 16 ]]
then
	if [[ -z "${FIX_REPMASK_USELAFIX_PATH}" ]]
	then 
		(>&2 echo "WARNING - Variable FIX_REPMASK_USELAFIX_PATH is not set.Try to use default path: patchedReads_dalign")
		FIX_REPMASK_USELAFIX_PATH="patchedReads_dalign"
	fi
		
	mkdir -p ${ASSMEBLY_DIR}/${FIX_REPMASK_USELAFIX_PATH}
	cd ${ASSMEBLY_DIR}/${FIX_REPMASK_USELAFIX_PATH}
	${SUBMIT_SCRIPTS_PATH}/createAndSubmitMarvelSlurmJobs.sh ${realPathConfigFile} ${currentPhase} ${currentStep} ${Id}
	cd ${cwd}
fi
