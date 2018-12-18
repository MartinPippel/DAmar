#!/bin/bash -e

configFile=$1
currentStep=$2
slurmID=$3

if [[ ! -f ${configFile} ]]
then 
    (>&2 echo "cannot access config file ${configFile}")
    exit 1
fi

source ${configFile}

gsize=${GSIZE}
i=$((${#GSIZE}-1))
if [[ "${GSIZE: -1}" =~ [gG] ]]
then
 gsize=$((${GSIZE:0:$i}*1000*1000*1000))
fi
if [[ "${GSIZE: -1}" =~ [mM] ]]
then
 gsize=$((${GSIZE:0:$i}*1000*1000))
fi
if [[ "${GSIZE: -1}" =~ [kK] ]]
then
 gsize=$((${GSIZE:0:$i}*1000))
fi

if [[ -z ${BIONANO_PATH} || ! -d ${BIONANO_PATH} ]]
then
	(>&2 echo "Variable BIONANO_PATH must be set to a proper bionano solve installation directory!!")
    exit 1
fi

function setBionanoOptions()
{
	if [[ $1 -eq 1 ]] ## single enzyme workflow
	then
		BIONANO_OPT="-f"
		
		if [[ -z ${SC_BIONANO_CONFLICTLEVEL_GENOMEMAPS} ]]
		then
			(>&2 echo "[WARNING] - Variable SC_BIONANO_CONFLICTLEVEL_GENOMEMAPS not set. Use 2, i.e. cut cmap at conflict!!")
			SC_BIONANO_CONFLICTLEVEL_GENOMEMAPS=2
		fi
		
		if [[ -z ${SC_BIONANO_CONFLICTLEVEL_SEQUENCE} ]]
		then
			(>&2 echo "[WARNING] - Variable SC_BIONANO_CONFLICTLEVEL_SEQUENCE not set. Use 2, i.e. cut contig at conflict!!")
			SC_BIONANO_CONFLICTLEVEL_SEQUENCE=2
		fi
		
		BIONANO_OPT="${BIONANO_OPT} -B ${SC_BIONANO_CONFLICTLEVEL_GENOMEMAPS} -N ${SC_BIONANO_CONFLICTLEVEL_SEQUENCE}"
		
		if [[ -z ${SC_BIONANO_MOLCULES_1} || -z ${SC_BIONANO_ASSEMBLY_1} || -z ${SC_BIONANO_ASSEMBLYSCRIPT_1} || -z ${SC_BIONANO_ASSEMBLYOPTARGS_1} || -z ${SC_BIONANO_ASSEMBLY_NOISE1} ]]
		then 
			(>&2 echo "[WARNING] - To align molecules to hybrid scaffolds (-x) and to generate a chimeric quality score all of the following variables have to be set properly:")
			(>&2 echo "            SC_BIONANO_MOLCULES_1, SC_BIONANO_ASSEMBLY_1, SC_BIONANO_ASSEMBLYSCRIPT_1, SC_BIONANO_ASSEMBLYOPTARGS_1, SC_BIONANO_ASSEMBLY_NOISE1"
		else if [[ ! -f "${SC_BIONANO_MOLCULES_1}" || ! -f "${SC_BIONANO_ASSEMBLY_1}" || ! -f "${SC_BIONANO_ASSEMBLYSCRIPT_1}" || ! -f "${SC_BIONANO_ASSEMBLYOPTARGS_1}" || ! -f "${SC_BIONANO_ASSEMBLY_NOISE1}" ]]
		then 
			(>&2 echo "[WARNING] - To align molecules to hybrid scaffolds (-x) and to generate a chimeric quality score all of the following variables have to be set properly:")
			(>&2 echo "            SC_BIONANO_MOLCULES_1, SC_BIONANO_ASSEMBLY_1, SC_BIONANO_ASSEMBLYSCRIPT_1, SC_BIONANO_ASSEMBLYOPTARGS_1, SC_BIONANO_ASSEMBLY_NOISE1"
		else
			BIONANO_OPT="${BIONANO_OPT} -x -y"
			BIONANO_OPT="${BIONANO_OPT} -m ${SC_BIONANO_MOLCULES_1}"
			BIONANO_OPT="${BIONANO_OPT} -p ${SC_BIONANO_ASSEMBLYSCRIPT_1}"
			BIONANO_OPT="${BIONANO_OPT} -q ${SC_BIONANO_ASSEMBLYOPTARGS_1}"
			BIONANO_OPT="${BIONANO_OPT} -e ${SC_BIONANO_ASSEMBLY_NOISE_1}"
		fi
	fi
}

myTypes=("01_BNscaffold")
if [[ ${SC_BIONANO_TYPE} -eq 0 ]]
then 
    ### 01_BNscaffold
    if [[ ${currentStep} -eq 1 ]]
    then
        ### clean up plans 
        for x in $(ls bionano_01_*_*_${CONT_DB}.${slurmID}.* 2> /dev/null)
        do            
            rm $x
        done
        
        # check if cmap(s) exist
        if [[ -z "${SC_BIONANO_ASSEMBLY_1}" || ! -f ${SC_BIONANO_ASSEMBLY_1} ]]
        then
        	(>&2 echo "ERROR - set SC_BIONANO_ASSEMBLY_1 to proper bionano de novo assembly cmap file")
        	exit 1
   		fi
   		
   		if [[ -z ${SC_BIONANO_ENZYME_1} || ! ${SC_BIONANO_ENZYME_1} =~ ^(DLE1|BSSSI|BSPQI)$ ]]
   		then 
   			(>&2 echo "ERROR - set SC_BIONANO_ENZYME_1 to proper bionano enzyme type [DLE1, BSPQI, BSSSI]")
        	exit 1	
   		fi
   	
   		if [[ -n "${SC_BIONANO_ASSEMBLY_2}" && -f ${SC_BIONANO_ASSEMBLY_2} && ${SC_BIONANO_ASSEMBLY_1} != ${SC_BIONANO_ASSEMBLY_2} ]]
        then
        	TWO_ENZYME_WORKFLOW=1
        	if [[ -z ${SC_BIONANO_ENZYME_2} || ! ${SC_BIONANO_ENZYME_2} =~ ^(DLE1|BSSSI|BSPQI)$ || ${SC_BIONANO_ENZYME_1} == ${SC_BIONANO_ENZYME_2} ]]
   			then 
   				(>&2 echo "ERROR - set SC_BIONANO_ENZYME_2 to proper bionano enzyme type [DLE1, BSPQI, BSSSI]")
   				(>&2 echo "        and SC_BIONANO_ENZYME_1 must be differen from SC_BIONANO_ENZYME_2")
        		exit 1	
   			fi        	
    	fi
    	       		   		
   		if [[ -n ${SC_BIONANO_OUTDIR} ]]
   		then 
   			(>&2 echo "ERROR - Variable SC_BIONANO_OUTDIR must be set!")
	        exit 1	
   		fi
   		
   		if [[ -n ${SC_BIONANO_RUNID} ]]
   		then 
   			SC_BIONANO_RUNID=1
   		fi 
   		
   		if [[ ! -f ${SC_BIONANO_REF} ]]
   		then
   			(>&2 echo "ERROR - set SC_BIONANO_REF to proper reference fasta file")
        	exit 1	
   		fi
   		
   		echo "if [[ -d ${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID} ]]; then mv ${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID} ${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}_$(date '+%Y-%m-%d_%H-%M-%S'); fi && mkdir -p ${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}" > bionano_01_BNscaffold_single_${CONT_DB}.${slurmID}.plan
   		echo "mkdir -p ${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/ref" >> bionano_01_BNscaffold_single_${CONT_DB}.${slurmID}.plan
   		echo "mkdir -p ${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/cmaps" >> bionano_01_BNscaffold_single_${CONT_DB}.${slurmID}.plan
   		echo "mkdir -p ${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/config" >> bionano_01_BNscaffold_single_${CONT_DB}.${slurmID}.plan
   		echo "mkdir -p ${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/out" >> bionano_01_BNscaffold_single_${CONT_DB}.${slurmID}.plan
   		echo "ln -s -r ${SC_BIONANO_REF} ${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/ref" >> bionano_01_BNscaffold_single_${CONT_DB}.${slurmID}.plan
   		echo "ln -s -r ${SC_BIONANO_ASSEMBLY_1} ${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/cmaps" >> bionano_01_BNscaffold_single_${CONT_DB}.${slurmID}.plan
   	
   		if [[ ${TWO_ENZYME_WORKFLOW} ]]
		then
			echo "ln -s -r ${SC_BIONANO_ASSEMBLY_2} ${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/cmaps" >> bionano_01_BNscaffold_single_${CONT_DB}.${slurmID}.plan
   		
   		
   		else
   			# set path hybrid scaffolder 
   			HYBSCAF=$(ls ${BIONANO_PATH}/HybridScaffold/????????/hybridScaffold.pl)
   			if [[ ! -f ${HYBSCAF} ]]
   			then 
   				(>&2 echo "ERROR - Cannot find hybridScaffold.pl. Should be here: ${BIONANO_PATH}/HybridScaffold/????????/hybridScaffold.pl")
        		exit 1
   			fi
   			HYBSCAF_PATH=$(dirname ${HYBSCAF})
   			# set reference variable
   			REF="${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/ref/$(basename ${SC_BIONANO_REF})"
   			# set cmap variable
   			CMAP=${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/cmaps/$(basename ${SC_BIONANO_ASSEMBLY_1})
   			# set hybrid config file
   			if [[ ${SC_BIONANO_ENZYME_1} == "DLE1" ]] 
   			then 
   				echo "cp ${HYBSCAF_PATH}/hybridScaffold_DLE1_config.xml ${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/config/"
   				HYB_CONF=${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/config/hybridScaffold_DLE1_config.xml
   			elif [[ ${SC_BIONANO_ENZYME_1} == "BSPQI" ]]
   			then
   				echo "sed -e "s:CTTAAG:BspQI:" ${HYBSCAF_PATH}/hybridScaffold_config.xml ${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/config/hybridScaffold_BSPQI_config.xml"
				HYB_CONF=${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/config/hybridScaffold_BSPQI_config.xml
   			elif [[ ${SC_BIONANO_ENZYME_1} == "BSSSI" ]]
   			then
   				echo "sed -e "s:CTTAAG:BssSI:" ${HYBSCAF_PATH}/hybridScaffold_config.xml ${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/config/hybridScaffold_BSSSI_config.xml"
				HYB_CONF=${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/config/hybridScaffold_BSSSI_config.xml				   				
   			fi
   			# set refaligner path 
   			REFALN=$(find ${BIONANO_PATH}/RefAligner/*/avx -name "RefAligner")
   			if [[ $(echo -e "${REFALN}" | wc -l ) -ne 1 ]]
   			then
   				(>&2 echo "ERROR - Cannot find RefAligner binary. Should be here: ${BIONANO_PATH}/RefAligner/*/avx/RefAligner")
        		exit 1
   			fi 
   			# set output directory
   			OUT="${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/out"
   			
   			setBionanoOptions 1
   			
   			echo "perl ${HYBSCAF} -n ${REF} -b ${CMAP} -c ${HYB_CONF} -r ${REFALN} -o ${OUT} -f -B 2 -N 2 -x -y ${BIONANO_OPT}" > bionano_01_BNscaffold_single_${CONT_DB}.${slurmID}.plan
   			echo "${REFALN}" > bionano_01_BNscaffold_single_${CONT_DB}.${slurmID}.version
   			echo "$(perl ${HYBSCAF} -v)" >> bionano_01_BNscaffold_single_${CONT_DB}.${slurmID}.version
   		fi
    else
        (>&2 echo "step ${currentStep} in SC_BIONANO_TYPE ${SC_BIONANO_TYPE} not supported")
        (>&2 echo "valid steps are: ${myTypes[${SC_BIONANO_TYPE}]}")
        exit 1            
    fi    		
else
    (>&2 echo "unknown SC_BIONANO_TYPE ${SC_BIONANO_TYPE}")
    (>&2 echo "supported types")
    x=0; while [ $x -lt ${#myTypes[*]} ]; do (>&2 echo "${myTypes[${x}]}"); done 
    exit 1
fi

exit 0