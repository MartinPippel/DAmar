#!/bin/bash 

config=$1
phase=$2

if [[ ! -f ${config} ]]
then 
  (>&2 echo "config ${config} not available")
  exit 1
fi

source ${config}

if [[ -z "${PROJECT_ID}" ]]
then 
    (>&2 echo "ERROR - You have to specify a project id. Set variable PROJECT_ID")
    exit 1
fi

if [[ -z ${QUAST_PATH} || ! -f ${QUAST_PATH}/quast.py ]]
then
	(>&2 echo "Variable QUAST_PATH must be set to a proper quast installation directory!!")
    exit 1
fi

if [[ ! -f ${SEQKIT_PATH} ]]
then
	(>&2 echo "Variable SEQKIT_PATH must be set to a proper seqkit binary!!")
    exit 1
fi

mkdir -p stats

#### create fasta stats (haploid/diploid size, ng50s, etc)

db=${DB%.db}
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

if [[ ${phase} -eq 6 ]] ## raw assembly stats  (last step in touring)
then 
	if [[ -d ${FIX_FILT_OUTDIR}/tour ]]
	then
		p=stats/contigs/${FIX_FILT_OUTDIR}/raw
		if [[ -d ${p} ]]
		then
			mv ${p} ${p}_$(date '+%Y-%m-%d_%H-%M-%S')	
		fi
		mkdir -p ${p}
		
		if [[ -f ${p}/coverage.hist ]]
		then 
			rm ${p}/coverage.hist
		fi
		# create effective coverage histogram - DAScover(based on RAW_DB) 
		if [[ -f ../../${COVERAGE_DIR}/effectiveCov_${RAW_DAZZ_DB%.db}_1_forBlock${RAW_DASCOVER_DALIGNER_FORBLOCK}.txt ]]
		then 
			echo "DAScover $(grep -e "Coverage is estimated at" ../../${COVERAGE_DIR}/effectiveCov_${RAW_DAZZ_DB%.db}_1_forBlock${RAW_DASCOVER_DALIGNER_FORBLOCK}.txt | awk '{print $NF}')" > ${p}/coverage.hist
			grep -e "%" ../../${COVERAGE_DIR}/effectiveCov_${RAW_DAZZ_DB%.db}_1_forBlock${RAW_DASCOVER_DALIGNER_FORBLOCK}.txt  | tr -d ":" | awk '{print $1" "$2}' | tac > ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_DAScover_b${RAW_DASCOVER_DALIGNER_FORBLOCK}.hist	
		fi
		
		# create effective coverage histogram - LArepeat(based on FIX_DB)
		for y in log_scrub_LArepeat_${RAW_DB}/*_1.out
		do
			if [[ -f ${y} ]]
			then
				grep -e "^COV " ${y%1.out}* | sort -k2,2n | awk -v c=0 -v cumCov=0 '{if ($2==c) {cumCov+=$4} else {print cumCov" "c; c=$2; cumCov=$4}} END {print cumCov" "c}' > ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_LArepeat_$(basename ${y%_1.out}).hist
				### tail -n +3 --- ignore counts at coverage 0 and 1
				echo "LArepeat_$(basename ${y%_1.out}) $(cat ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_LArepeat_$(basename ${y%_1.out}).hist | tail -n +3 | sort -k2,2n | tail -n1 | awk '{print $1}')" >> ${p}/coverage.hist	
			fi			
		done
		
		for y in ${FIX_FILT_OUTDIR}/tour/*.fasta
		do
			name=$(basename ${y%.fasta})
			sed -e "s:>path_:>${name}_:" $y  
		done > ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_r.fasta
		${SUBMIT_SCRIPTS_PATH}/splitDiploidAssembly.py ${PROJECT_ID}_${FIX_FILT_OUTDIR}_r ${gsize} ${p} ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_r.fasta  

		## create statistic files
		cat ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_r.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_r.stats
		cat ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_r.p.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_r.p.stats
		cat ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_r.a.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_r.a.stats
		
		## create header files
		grep -e ">" ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_r.fasta | awk '{print $1}' | sed -e 's:^>::' > ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_r.header
		grep -e ">" ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_r.p.fasta | awk '{print $1}' | sed -e 's:^>::' > ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_r.p.header
		grep -e ">" ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_r.a.fasta | awk '{print $1}' | sed -e 's:^>::' > ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_r.a.header
		
		## copy config file
		cp $config ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_$(date '+%Y-%m-%d_%H-%M-%S').config.sh
	else
		(>&2 echo "ERROR - directory ${FIX_FILT_OUTDIR}/tour not available")
  		exit 1
	fi	
elif [[ ${phase} -eq 7 ]] ## marvel corrected assembly stats  (last step in correction)
then
	if [[ -d ${FIX_FILT_OUTDIR}/${COR_DIR}/contigs ]]
	then
		p=stats/contigs/${FIX_FILT_OUTDIR}/corr
		if [[ -d ${p} ]]
		then
			mv ${p} ${p}_$(date '+%Y-%m-%d_%H-%M-%S')	
		fi
		mkdir -p ${p}
		
		for y in ${FIX_FILT_OUTDIR}/${COR_DIR}/contigs/*.fasta
		do
			cat $y  
		done > ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_c.fasta
		${SUBMIT_SCRIPTS_PATH}/splitDiploidAssembly.py ${PROJECT_ID}_${FIX_FILT_OUTDIR}_c ${gsize} ${p} ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_c.fasta  

		## create statistic files
		cat ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_c.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_c.stats
		cat ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_c.p.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_c.p.stats
		cat ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_c.a.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_c.a.stats
		
		## copy config file
		cp $config ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_$(date '+%Y-%m-%d_%H-%M-%S').config.sh
	else
		(>&2 echo "ERROR - directory ${FIX_FILT_OUTDIR}/${COR_DIR}/contigs not available")
  		exit 1
	fi
	
elif [[ ${phase} -eq 8 ]] ## marvel corrected assembly stats after contig analysis  (last step after CTanalysis)
then
	if [[ -d ${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${COR_CONTIG_CTANALYZE_DIR}/classified ]]
	then
		p=stats/contigs/${FIX_FILT_OUTDIR}/haploSplit
		if [[ -d ${p} ]]
		then
			mv ${p} ${p}_$(date '+%Y-%m-%d_%H-%M-%S')	
		fi
		mkdir -p ${p}
		mkdir -p ${p}/forArrow
		
		prim=${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${COR_CONTIG_CTANALYZE_DIR}/classified/asm.p.fa
		alt=${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${COR_CONTIG_CTANALYZE_DIR}/classified/asm.a.fa
		crap=${FIX_FILT_OUTDIR}/${ANALYZE_DIR}/${COR_CONTIG_CTANALYZE_DIR}/classified/asm.c.fa
	
		if [[ ! -f ${prim} ]]
		then
			(>&2 echo "ERROR - primary contig file missing: ${prim}")
  			exit 1  			
  		elif [[ ! -f ${alt} ]]
		then
			(>&2 echo "ERROR - alt contig file missing: ${alt}")
  			exit 1
		elif [[ ! -f ${crap} ]]
		then
			(>&2 echo "ERROR - crap contig file missing: ${crap}")
  			exit 1  	
		fi
		
		cp ${prim} ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_h.p.fasta
		cp ${alt} ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_h.a.fasta
		cp ${crap} ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_h.c.fasta
		
		cat ${prim} ${alt} > ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_h.fasta
				
		grep -e ">" ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_h.fasta | awk '{print $1}' | sed -e 's:^>::' > ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_h.header
		grep -e ">" ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_h.p.fasta | awk '{print $1}' | sed -e 's:^>::' > ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_h.p.header
		grep -e ">" ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_h.a.fasta | awk '{print $1}' | sed -e 's:^>::' > ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_h.a.header				
		grep -e ">" ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_h.c.fasta | awk '{print $1}' | sed -e 's:^>::' > ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_h.c.header
	
		ln -s -r ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_h.fasta ${p}/forArrow/${PROJECT_ID}_${FIX_FILT_OUTDIR}_h.fasta
		
		ln -s -r -f ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_h.p.header ${p}/forArrow/${PROJECT_ID}_${FIX_FILT_OUTDIR}_h.p.header
		ln -s -r -f ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_h.a.header ${p}/forArrow/${PROJECT_ID}_${FIX_FILT_OUTDIR}_h.a.header
	
		## create statistic files
		cat ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_h.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_h.stats
		cat ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_h.p.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_h.p.stats
		cat ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_h.a.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_h.a.stats
		cat ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_h.c.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_h.c.stats
		
		## copy config file
		cp $config ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_$(date '+%Y-%m-%d_%H-%M-%S').config.sh
		
		## copy assembly relationship stats as well as properties
		cp ${prim%.fa}.stats ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_h.p.cstats
		cp ${alt%.fa}.stats ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_h.a.cstats
		cp ${crap%.fa}.stats ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_h.c.cstats
		
		## FOR TESTING PURPOSE ONLY
		mkdir -p ${p}/filter
		## add some primary contig filter based on repeats
		header="count cumBases $(head -n 1 ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_h.p.cstats)"  
		for x in $(seq 65 1 80)
		do
			echo "${header}" > ${p}/filter/${PROJECT_ID}_${FIX_FILT_OUTDIR}_h.p.excludeP${x}RepeatContigs.cstats
			tail -n +2 ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_h.p.cstats | awk -v t=${x} '{if ($3>t || $4 > t) {s+=$2; c+=1; print c"\t"s"\t"$0}}' >> ${p}/filter/${PROJECT_ID}_${FIX_FILT_OUTDIR}_h.p.excludeP${x}RepeatContigs.cstats
			tail -n +2 ${p}/filter/${PROJECT_ID}_${FIX_FILT_OUTDIR}_h.p.excludeP${x}RepeatContigs.cstats | awk '{print $8"_"$7}' > ${p}/filter/${PROJECT_ID}_${FIX_FILT_OUTDIR}_h.p.excludeP${x}RepeatContigs.clist 
			${SEQKIT_PATH} grep -f ${p}/filter/${PROJECT_ID}_${FIX_FILT_OUTDIR}_h.p.excludeP${x}RepeatContigs.clist ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_h.p.fasta > ${p}/filter/${PROJECT_ID}_${FIX_FILT_OUTDIR}_h.pEP${x}.fasta
			${SEQKIT_PATH} grep -v -f ${p}/filter/${PROJECT_ID}_${FIX_FILT_OUTDIR}_h.p.excludeP${x}RepeatContigs.clist ${p}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_h.p.fasta > ${p}/filter/${PROJECT_ID}_${FIX_FILT_OUTDIR}_h.pIP${x}.fasta
			cat ${p}/filter/${PROJECT_ID}_${FIX_FILT_OUTDIR}_h.pIP${x}.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${p}/filter/${PROJECT_ID}_${FIX_FILT_OUTDIR}_h.pIP${x}.stats
			cat ${p}/filter/${PROJECT_ID}_${FIX_FILT_OUTDIR}_h.pEP${x}.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${p}/filter/${PROJECT_ID}_${FIX_FILT_OUTDIR}_h.pEP${x}.stats 
		done
	else
		(>&2 echo "ERROR - directory ${FIX_FILT_OUTDIR}/correction/contigs not available")
  		exit 1
	fi
elif [[ ${phase} -eq 9 ]] ## assembly stats after PacBio Arrow Correction 
then
	if [[ -d ${FIX_FILT_OUTDIR}/arrow_${PB_ARROW_RUNID} ]]
	then
		
		arrowPath=stats/contigs/${FIX_FILT_OUTDIR}/arrow_${PB_ARROW_RUNID}

		fext="a"
		prevExt=$(basename ${PB_ARROW_REFFASTA%.fasta} | awk -F '[_.]' '{print $(NF-1)}')	
		ext=$(basename ${PB_ARROW_REFFASTA%.fasta} | awk -F '[_.]' '{print $NF}')	
		
		if [[ ! -f ${PB_ARROW_REFFASTA} ]]
        then 
	     (>&2 echo "ERROR - Reference file name not set PB_ARROW_P_HEADER=${PB_ARROW_REFFASTA}.")
    	    exit 1        			
        fi
        
        if [[ -d ${arrowPath} ]]
		then
			mv ${arrowPath} ${arrowPath}_$(date '+%Y-%m-%d_%H-%M-%S')	
		fi
        mkdir -p ${arrowPath}        
		mkdir -p ${arrowPath}/forArrow
	
		## primary 
		if [[ -n ${PB_ARROW_P_HEADER} && -f ${PB_ARROW_P_HEADER} ]]
		then
			for z in $(cat ${PB_ARROW_P_HEADER})
			do
				name=$(echo ${z} | awk -F \_ '{$NF=""; OFS="_"; print $0}')
				pathID=$(echo ${z} | awk -F \_ '{print $NF}')				
				
				if [[ ! -d ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${name}${pathID} ]]
				then 
					echo "WARNING - Missing directory ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${name}${pathID}." >> ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.p.elog
					continue
				fi  
				
				if [[ ! -f ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${name}${pathID}/ALL_${name}${pathID}.arrow.fa ]]
				then 
				 	echo "WARNING - Arrow failed. Missing file ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${name}${pathID}/ALL_${name}${pathID}.arrow.fa." >> ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.p.elog
					continue
				fi  
				
				cat ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${name}${pathID}/ALL_${name}${pathID}.arrow.fa	| sed -e "s:[|_]arrow::g"					        		
			done > ${arrowPath}/${xPROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.p.fasta
			${QUAST_PATH}/quast.py -t 1 -s -e --fast --est-ref-size ${gsize} -o ${arrowPath}/${xPROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.p ${arrowPath}/${xPROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.p.fasta
			if [[ -n ${PB_ARROW_BGZIP} && ${PB_ARROW_BGZIP} -gt 1 ]]
			then
				${CONDA_BASE_ENV} && bgzip -@${PB_ARROW_BGZIP} -c ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.p.fasta > ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.p.fa.gz && conda deactivate
			else 
				gzip -c ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.p.fasta > ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.p.fa.gz
			fi
		fi
		if [[-n ${PB_ARROW_P_HEADER} && -f ${PB_ARROW_P_HEADER} ]]
		then 
			## alternate
			for z in $(cat ${PB_ARROW_P_HEADER})
			do
				name=$(echo ${z} | awk -F \_ '{$NF=""; OFS="_"; print $0}')
				pathID=$(echo ${z} | awk -F \_ '{print $NF}')
				
				if [[ ! -d ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${name}${pathID} ]]
				then 
					echo "WARNING - Missing directory ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${name}${pathID}." >> ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.a.elog
					continue
				fi  
				
				if [[ ! -f ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${name}${pathID}/ALL_${name}${pathID}.arrow.fa ]]
				then 
				 	echo "WARNING - Arrow failed. Missing file ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${name}${pathID}/ALL_${name}${pathID}.arrow.fa." >> ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.a.elog
					continue
				fi  
				
				cat ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${name}${pathID}/ALL_${name}${pathID}.arrow.fa	| sed -e "s:[|_]arrow::g"					        		
			done > ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.a.fasta
			${QUAST_PATH}/quast.py -t 1 -s -e --fast --est-ref-size ${gsize} -o ${arrowPath}/${xPROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.a ${arrowPath}/${xPROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.a.fasta
			if [[ -n ${PB_ARROW_BGZIP} && ${PB_ARROW_BGZIP} -gt 1 ]]
			then
				${CONDA_BASE_ENV} && bgzip -@${PB_ARROW_BGZIP} -c ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.a.fasta > ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.a.fa.gz && conda deactivate
			else 
				gzip -c ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.a.fasta > ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.a.fa.gz
			fi
		fi
		## all
		for z in $(grep -e ">" ${PB_ARROW_REFFASTA} | awk -F '[> ]' '{print $2}')
		do
			name=$(echo ${z} | awk -F \_ '{$NF=""; OFS="_"; print $0}')
			pathID=$(echo ${z} | awk -F \_ '{print $NF}')
			
			if [[ ! -d ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${name}${pathID} ]]
			then 
				echo "WARNING - Missing directory ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${name}${pathID}." >> ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.${ext}.elog
				continue
			fi  
			
			if [[ ! -f ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${name}${pathID}/ALL_${name}${pathID}.arrow.fa ]]
			then 
			 	echo "WARNING - Arrow failed. Missing file ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${name}${pathID}/ALL_${name}${pathID}.arrow.fa." >> ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.${ext}.elog
				continue
			fi
			cat ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${name}${pathID}/ALL_${name}${pathID}.arrow.fa	| sed -e "s:[|_]arrow::g"
		done > ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.${ext}.fasta
		${QUAST_PATH}/quast.py -t 1 -s -e --fast --est-ref-size ${gsize} -o ${arrowPath}/${xPROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext} ${arrowPath}/${xPROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.fasta
		if [[ -n ${PB_ARROW_BGZIP} && ${PB_ARROW_BGZIP} -gt 1 ]]
		then
			${CONDA_BASE_ENV} && bgzip -@${PB_ARROW_BGZIP} -c ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.${ext}.fasta > ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.${ext}.fa.gz && conda deactivate
		else
			gzip -c ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.${ext}.fasta > ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.${ext}.fa.gz
       	fi
       
        if [[ -s ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.${ext}.fasta ]]
        then 
        	cat ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.${ext}.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.${ext}.stats
        	allBases=$(cat ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.${ext}.fasta | awk '{ if ($1 ~ /^>/) ; else print $0;}' | tr -d "\n" | wc -m)
        	arrowedBases=$(cat ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.${ext}.fasta | awk '{ if ($1 ~ /^>/) ; else print $0;}' | tr -d "\nacgt" | wc -m)
        	frac=$(echo "scale=4;${arrowedBases}*100/${allBases}" | bc)
        	echo "allBases ${allBases} arrowedBases ${arrowedBases} arrowedFraction ${frac}%" >> ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.${ext}.stats
        	${SUBMIT_SCRIPTS_PATH}/trimLowerCaseTips.py ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.${ext}.fasta ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.${ext}T.fasta ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.${ext}T.failed.fasta
        	cat ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.${ext}T.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.${ext}T.stats
			allBases=$(cat ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.${ext}T.fasta | awk '{ if ($1 ~ /^>/) ; else print $0;}' | tr -d "\n" | wc -m)
        	arrowedBases=$(cat ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.${ext}T.fasta | awk '{ if ($1 ~ /^>/) ; else print $0;}' | tr -d "\nacgt" | wc -m)
        	frac=$(echo "scale=4;${arrowedBases}*100/${allBases}" | bc)
        	echo "allBases ${allBases} arrowedBases ${arrowedBases} arrowedFraction ${frac}%" >> ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.${ext}T.stats
    	fi
        	
        if [[ -s ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.p.fasta ]]
        then
        	cat ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.p.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.p.stats
        	allBases=$(cat ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.p.fasta | awk '{ if ($1 ~ /^>/) ; else print $0;}' | tr -d "\n" | wc -m)
        	arrowedBases=$(cat ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.p.fasta | awk '{ if ($1 ~ /^>/) ; else print $0;}' | tr -d "\nacgt" | wc -m)
        	frac=$(echo "scale=4;${arrowedBases}*100/${allBases}" | bc)
        	echo "allBases ${allBases} arrowedBases ${arrowedBases} arrowedFraction ${frac}%" >> ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.p.stats
        	${SUBMIT_SCRIPTS_PATH}/trimLowerCaseTips.py ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.p.fasta ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.pT.fasta ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.pT.failed.fasta
        	cat ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.pT.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.pT.stats
			allBases=$(cat ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.pT.fasta | awk '{ if ($1 ~ /^>/) ; else print $0;}' | tr -d "\n" | wc -m)
        	arrowedBases=$(cat ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.pT.fasta | awk '{ if ($1 ~ /^>/) ; else print $0;}' | tr -d "\nacgt" | wc -m)
        	frac=$(echo "scale=4;${arrowedBases}*100/${allBases}" | bc)
        	echo "allBases ${allBases} arrowedBases ${arrowedBases} arrowedFraction ${frac}%" >> ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.pT.stats
    	fi 
       	
        if [[ -s ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.a.fasta ]]
        then
        	cat ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.a.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.a.stats
        	allBases=$(cat ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.a.fasta | awk '{ if ($1 ~ /^>/) ; else print $0;}' | tr -d "\n" | wc -m)
        	arrowedBases=$(cat ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.a.fasta | awk '{ if ($1 ~ /^>/) ; else print $0;}' | tr -d "\nacgt" | wc -m)
        	frac=$(echo "scale=4;${arrowedBases}*100/${allBases}" | bc )
        	echo "allBases ${allBases} arrowedBases ${arrowedBases} arrowedFraction ${frac}%" >> ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.a.stats
        	${SUBMIT_SCRIPTS_PATH}/trimLowerCaseTips.py ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.a.fasta ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.aT.fasta ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.aT.failed.fasta
        	cat ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.aT.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.aT.stats
			allBases=$(cat ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.aT.fasta | awk '{ if ($1 ~ /^>/) ; else print $0;}' | tr -d "\n" | wc -m)
        	arrowedBases=$(cat ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.aT.fasta | awk '{ if ($1 ~ /^>/) ; else print $0;}' | tr -d "\nacgt" | wc -m)
        	frac=$(echo "scale=4;${arrowedBases}*100/${allBases}" | bc)        	
        	echo "allBases ${allBases} arrowedBases ${arrowedBases} arrowedFraction ${frac}%" >> ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.aT.stats
    	fi	 
    	
    	## prepare files for another arrow run
    	if [[ -s ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.p.fasta ]]
    	then
    		grep -e ">" ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.p.fasta | awk '{print $1}' | sed -e 's:^>::' > ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.p.header
    		ln -s -r -f ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.p.header ${arrowPath}/forArrow/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.p.header
    	fi
    	
    	if [[ -s ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${prevExt}${fext}.a.fasta ]]
    	then
			grep -e ">" ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.a.fasta | awk '{print $1}' | sed -e 's:^>::' > ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.a.header
			ln -s -r -f ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.a.header ${arrowPath}/forArrow/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.a.header
		fi				
		
		ln -s -r -f ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.fasta ${arrowPath}/forArrow/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${prevExt}${fext}.fasta
	else
		(>&2 echo "ERROR - directory ${FIX_FILT_OUTDIR}/arrow_${PB_ARROW_RUNID} not available")
  		exit 1
	fi
elif [[ ${phase} -eq 11 ]] ## after freebayes polishing 
then
	if [[ -d ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID} ]]
	then
		freebayesPath=stats/contigs/${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}
		fext="f"
		prevExt=$(basename ${CT_FREEBAYES_REFFASTA%.fasta} | awk -F '[_.]' '{print $(NF-1)}')		
		
		if [[ -d ${freebayesPath} ]]
		then
			mv ${freebayesPath} ${freebayesPath}_$(date '+%Y-%m-%d_%H-%M-%S')	
		fi
		mkdir -p ${freebayesPath}
		
		fname="${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/${PROJECT_ID}_10x.fasta"
		
		if [[ ! -f ${fname} ]]
		then
			(>&2 echo "ERROR - Cannot find polished file ${fname}")
  			exit 1
		fi
		
		cat ${fname} > ${freebayesPath}/${PROJECT_ID}_${CT_FREEBAYES_OUTDIR}_${prevExt}${fext}.p.fasta
		cp ${fname%.fasta}.numvar ${freebayesPath}/${PROJECT_ID}_${CT_FREEBAYES_OUTDIR}_${prevExt}${fext}.p.numvar
		gzip -c ${freebayesPath}/${PROJECT_ID}_${CT_FREEBAYES_OUTDIR}_${prevExt}${fext}.p.fasta > ${freebayesPath}/${PROJECT_ID}_${CT_FREEBAYES_OUTDIR}_${prevExt}${fext}.p.fa.gz
		cat ${freebayesPath}/${PROJECT_ID}_${CT_FREEBAYES_OUTDIR}_${prevExt}${fext}.p.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${freebayesPath}/${PROJECT_ID}_${CT_FREEBAYES_OUTDIR}_${prevExt}${fext}.p.stats
		cp ${config} ${freebayesPath}/$(date '+%Y-%m-%d_%H-%M-%S')_$(basename ${config})
	else
		(>&2 echo "ERROR - directory ${FIX_FILT_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID} not available")
  		exit 1
	fi
elif [[ ${phase} -eq 12 ]] ## 10x scaffolding 
then
	if [[ ${SC_10X_TYPE} -eq 0 ]] && [[ -d ${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID} ]]
	then 
		
		scaff10xPath=stats/contigs/${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID}
		if [[ -d ${scaff10xPath} ]]
		then
			mv ${scaff10xPath} ${scaff10xPath}_$(date '+%Y-%m-%d_%H-%M-%S')	
		fi
		mkdir -p ${scaff10xPath}
		
		prevExt=$(basename ${SC_10X_REF%.fasta} | awk -F '[_.]' '{print $(NF-1)}')
        scaffdir="${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID}"
        # scaff10x step 2
        inputScaffold0="${scaffdir}/${PROJECT_ID}_${SC_10X_OUTDIR}_${prevExt}b.p.fasta"
        inputScaffold1="${scaffdir}/${PROJECT_ID}_${SC_10X_OUTDIR}_${prevExt}x.p.fasta" 
    	inputScaffold2="${scaffdir}/${PROJECT_ID}_${SC_10X_OUTDIR}_${prevExt}bx.p.fasta"
    	# scaff10x step 3a - break10x
    	inputScaffold3="${scaffdir}/${PROJECT_ID}_${SC_10X_OUTDIR}_${prevExt}xb.p.fasta" 
    	inputScaffold4="${scaffdir}/${PROJECT_ID}_${SC_10X_OUTDIR}_${prevExt}bxb.p.fasta"
    	# scaff10x step 3b - scaff10x
    	inputScaffold5="${scaffdir}/${PROJECT_ID}_${SC_10X_OUTDIR}_${prevExt}xx.p.fasta"
    	inputScaffold6="${scaffdir}/${PROJECT_ID}_${SC_10X_OUTDIR}_${prevExt}bxbx.p.fasta"
    	# scaff10x step 3c - break10x
    	inputScaffold7="${scaffdir}/${PROJECT_ID}_${SC_10X_OUTDIR}_${prevExt}xxb.p.fasta"
    	inputScaffold8="${scaffdir}/${PROJECT_ID}_${SC_10X_OUTDIR}_${prevExt}xbxb.p.fasta"
    	inputScaffold9="${scaffdir}/${PROJECT_ID}_${SC_10X_OUTDIR}_${prevExt}bxbxb.p.fasta"
    	
    	cp ${config} ${scaff10xPath}/$(date '+%Y-%m-%d_%H-%M-%S')_$(basename ${config})
    	for x in ${inputScaffold0} ${inputScaffold1} ${inputScaffold2} ${inputScaffold3} ${inputScaffold4} ${inputScaffold5} ${inputScaffold6} ${inputScaffold7} ${inputScaffold8} ${inputScaffold9}
    	do 
    		if [[ -f ${x} ]]
    		then
    			f=$(basename $x)
    			cp ${x} ${scaff10xPath}/
    			gzip -c ${scaff10xPath}/${f} > ${scaff10xPath}/${f%.fasta}.fa.gz
    			cat ${scaff10xPath}/${f} | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${scaff10xPath}/${f%.fasta}.stats
    			${QUAST_PATH}/quast.py -o ${scaff10xPath}/${f%.fasta} -t 1 -s -e --est-ref-size ${gsize} ${scaff10xPath}/${f}
    		else
    			(>&2 echo "WARNING assemblyStats.sh 12 - File ${x} is missing.")	
    		fi
    	done
    elif [[ ${SC_10X_TYPE} -eq 1 ]] && [[ -d ${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID} ]]
	then
		(>&2 echo "ERROR - 10x scaffolding type: ${SC_10X_TYPE}. Not implemented yet!")
	elif [[ ${SC_10X_TYPE} -eq 2 ]] && [[ -d ${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID} ]] ## the new scaff10x 4.0 pipeline 
	then
		scaff10xPath=stats/contigs/${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID}
		if [[ -d ${scaff10xPath} ]]
		then
			mv ${scaff10xPath} ${scaff10xPath}_$(date '+%Y-%m-%d_%H-%M-%S')	
		fi
		mkdir -p ${scaff10xPath}
		
		preName=$(basename ${SC_10X_REF%.fasta})
		prevExt=$(basename ${SC_10X_REF%.fasta} | awk -F '[_.]' '{print $(NF-1)}')
		prevSet=$(basename ${SC_10X_REF%.fasta} | awk -F '[_.]' '{print $(NF)}')
		
        scaffdir="${SC_10X_OUTDIR}/scaff10x_${SC_10X_RUNID}"
        # scaff10x step 2
        inputScaffold="${scaffdir}/${PROJECT_ID}_${SC_10X_OUTDIR}_${prevExt}x.${prevSet}.fasta"
		
		if [[ -f ${inputScaffold} ]]
		then
			f=$(basename $inputScaffold)
			cp ${inputScaffold} ${scaff10xPath}/
			cp ${inputScaffold}.agp ${scaff10xPath}/
			gzip -c ${scaff10xPath}/${f} > ${scaff10xPath}/${f%.fasta}.fa.gz
			cat ${scaff10xPath}/${f} | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${scaff10xPath}/${f%.fasta}.stats
			${QUAST_PATH}/quast.py -o ${scaff10xPath}/${f%.fasta} -t 1 -s -e --est-ref-size ${gsize} ${scaff10xPath}/${f}
			cp ${config} ${scaff10xPath}/$(date '+%Y-%m-%d_%H-%M-%S')_$(basename ${config})
		else
			(>&2 echo "WARNING assemblyStats.sh 12 - File ${inputScaffold} is missing.")	
			exit 1
		fi	
    elif [[ ${SC_10X_TYPE} -eq 3 ]] && [[ -d ${SC_10X_OUTDIR}/break10x_${SC_10X_RUNID} ]]
	then
		
		scaff10xPath=stats/contigs/${SC_10X_OUTDIR}/break10x_${SC_10X_RUNID}
		if [[ -d ${scaff10xPath} ]]
		then
			mv ${scaff10xPath} ${scaff10xPath}_$(date '+%Y-%m-%d_%H-%M-%S')	
		fi
		mkdir -p ${scaff10xPath}
		
		preName=$(basename ${SC_10X_REF%.fasta})
		prevExt=$(basename ${SC_10X_REF%.fasta} | awk -F '[_.]' '{print $(NF-1)}')
		prevSet=$(basename ${SC_10X_REF%.fasta} | awk -F '[_.]' '{print $(NF)}')
		
        scaffdir="${SC_10X_OUTDIR}/break10x_${SC_10X_RUNID}"
        # scaff10x step 2
        inputScaffold="${scaffdir}/${PROJECT_ID}_${SC_10X_OUTDIR}_${prevExt}b.${prevSet}.fasta"
		
		if [[ -f ${inputScaffold} ]]
		then
			f=$(basename $inputScaffold)
			cp ${inputScaffold} ${scaff10xPath}/
			cp ${inputScaffold%.fasta}.breaks ${scaff10xPath}/
			gzip -c ${scaff10xPath}/${f} > ${scaff10xPath}/${f%.fasta}.fa.gz
			cat ${scaff10xPath}/${f} | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${scaff10xPath}/${f%.fasta}.stats
			${QUAST_PATH}/quast.py -o ${scaff10xPath}/${f%.fasta} -t 1 -s -e --est-ref-size ${gsize} ${scaff10xPath}/${f}
			cp ${config} ${scaff10xPath}/$(date '+%Y-%m-%d_%H-%M-%S')_$(basename ${config})
		else
			(>&2 echo "WARNING assemblyStats.sh 12 - File ${inputScaffold} is missing.")	
			exit 1
		fi
	else
    (>&2 echo "ERROR - unknow 10x scaffolding type: ${SC_10X_TYPE} or directory ${SC_10X_OUTDIR}/[breal10x|scaff10x]s_${SC_10X_RUNID} not available")
		exit 1
	fi
elif [[ ${phase} -eq 13 ]] ## bionano scaffolding 
then
	if [[ -d ${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID} ]]
	then
		bionanoPath="stats/contigs/${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}"
		
		prevExt=$(basename ${SC_BIONANO_REF%.fasta} | awk -F '[_.]' '{print $(NF-1)}')
		cset=$(basename ${SC_BIONANO_REF%.fasta} | awk -F '[_.]' '{print $(NF)}')

		if [[ -d ${bionanoPath} ]]
		then
			mv ${bionanoPath} ${bionanoPath}_$(date '+%Y-%m-%d_%H-%M-%S')	
		fi				
		mkdir -p ${bionanoPath}
		
		PROJECT_ID_CAPS=$(echo ${PROJECT_ID} | awk '{print toupper($0)}')
		REF_NAME=$(basename ${SC_BIONANO_REF} | tr '.' '_')
		
		
		if [[ -f ${SC_BIONANO_REF_EXCLUSELIST} && -s ${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/exclude.fasta ]]
   		then
   			fext="E"
			cp  ${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/ref/exclude.fasta ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta
			cat ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.stats
			${QUAST_PATH}/quast.py -t 1 -s -e --fast --est-ref-size ${gsize} -o ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset} ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta
   		fi
		
		if [[ -n "${SC_BIONANO_ASSEMBLY_2}" && -f ${SC_BIONANO_ASSEMBLY_2} && ${SC_BIONANO_ASSEMBLY_1} != ${SC_BIONANO_ASSEMBLY_2} ]]
		then
			### two enzyme workflow 
			
			mkdir -p ${bionanoPath}/${SC_BIONANO_ENZYME_1}
			
			#--------------------------- enzyme-1
			fname="${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/out/${SC_BIONANO_ENZYME_1}/hybrid_scaffolds_M1/${PROJECT_ID_CAPS}_REFINEFINAL1_bppAdjust_cmap_${REF_NAME}_NGScontigs_HYBRID_SCAFFOLD_NCBI.fasta"
			agp="${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/out/${SC_BIONANO_ENZYME_1}/hybrid_scaffolds_M1/${PROJECT_ID_CAPS}_REFINEFINAL1_bppAdjust_cmap_${REF_NAME}_NGScontigs_HYBRID_SCAFFOLD.agp"
			if [[ ! -f ${fname} ]]
			then
				(>&2 echo "ERROR - Cannot find file ${fname}")
	  			exit 1
			fi
			# bionano scaffolds enzyme-1
			fext="n"
			cp ${fname} ${bionanoPath}/${SC_BIONANO_ENZYME_1}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta
			cp ${agp} ${bionanoPath}/${SC_BIONANO_ENZYME_1}
			gzip -c ${bionanoPath}/${SC_BIONANO_ENZYME_1}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta > ${bionanoPath}/${SC_BIONANO_ENZYME_1}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fa.gz
			cat ${bionanoPath}/${SC_BIONANO_ENZYME_1}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${bionanoPath}/${SC_BIONANO_ENZYME_1}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.stats
			${QUAST_PATH}/quast.py -t 1 -s -e --fast --est-ref-size ${gsize} -o ${bionanoPath}/${SC_BIONANO_ENZYME_1}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset} ${bionanoPath}/${SC_BIONANO_ENZYME_1}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta
			cp ${config} ${bionanoPath}/$(date '+%Y-%m-%d_%H-%M-%S')_$(basename ${config})
			cp ${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/out/${SC_BIONANO_ENZYME_1}/hybrid_scaffolds_M1/hybrid_scaffold_informatics_report.txt ${bionanoPath}/${SC_BIONANO_ENZYME_1}/
			# enzyme-1: cut bionano conflicts in input assembly
			fext="C" 
			cp ${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/out/${SC_BIONANO_ENZYME_1}/hybrid_scaffolds_M1/$(basename ${SC_BIONANO_REF}).cut.fasta ${bionanoPath}/${SC_BIONANO_ENZYME_1}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta
			cat ${bionanoPath}/${SC_BIONANO_ENZYME_1}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${bionanoPath}/${SC_BIONANO_ENZYME_1}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.stats
			cp ${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/out/${SC_BIONANO_ENZYME_1}/hybrid_scaffolds_M1/conflicts.txt ${bionanoPath}/${SC_BIONANO_ENZYME_1}/
			${QUAST_PATH}/quast.py -t 1 -e --fast --est-ref-size ${gsize} -o ${bionanoPath}/${SC_BIONANO_ENZYME_1}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset} ${bionanoPath}/${SC_BIONANO_ENZYME_1}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta
			# bionano not scaffolded 
			fname="${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/out/${SC_BIONANO_ENZYME_1}/hybrid_scaffolds_M1/${PROJECT_ID_CAPS}_REFINEFINAL1_bppAdjust_cmap_${REF_NAME}_NGScontigs_HYBRID_SCAFFOLD_NOT_SCAFFOLDED.fasta"
			if [[ ! -f ${fname} ]]
			then
				(>&2 echo "ERROR - Cannot find file ${fname}")
	  			exit 1
			fi
			fext="u" ### bionano not scaffolded
			cp ${fname} ${bionanoPath}/${SC_BIONANO_ENZYME_1}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta
			gzip -c ${bionanoPath}/${SC_BIONANO_ENZYME_1}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta > ${bionanoPath}/${SC_BIONANO_ENZYME_1}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fa.gz
			cat ${bionanoPath}/${SC_BIONANO_ENZYME_1}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${bionanoPath}/${SC_BIONANO_ENZYME_1}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.stats
			${QUAST_PATH}/quast.py -t 1 -s -e --fast --est-ref-size ${gsize} -o ${bionanoPath}/${SC_BIONANO_ENZYME_1}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset} ${bionanoPath}/${SC_BIONANO_ENZYME_1}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta
			# bionano combined (scaffolded + not scaffolded)
			fext="N" 
			cat ${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/out/${SC_BIONANO_ENZYME_1}/hybrid_scaffolds_M1/${PROJECT_ID_CAPS}_REFINEFINAL1_bppAdjust_cmap_${REF_NAME}_NGScontigs_HYBRID_SCAFFOLD_NCBI.fasta > ${bionanoPath}/${SC_BIONANO_ENZYME_1}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta
			cat ${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/out/${SC_BIONANO_ENZYME_1}/hybrid_scaffolds_M1/${PROJECT_ID_CAPS}_REFINEFINAL1_bppAdjust_cmap_${REF_NAME}_NGScontigs_HYBRID_SCAFFOLD_NOT_SCAFFOLDED.fasta >> ${bionanoPath}/${SC_BIONANO_ENZYME_1}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta
			gzip -c ${bionanoPath}/${SC_BIONANO_ENZYME_1}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta > ${bionanoPath}/${SC_BIONANO_ENZYME_1}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fa.gz
			cat ${bionanoPath}/${SC_BIONANO_ENZYME_1}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${bionanoPath}/${SC_BIONANO_ENZYME_1}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.stats
			${QUAST_PATH}/quast.py -t 1 -s -e --fast --est-ref-size ${gsize} -o ${bionanoPath}/${SC_BIONANO_ENZYME_1}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset} ${bionanoPath}/${SC_BIONANO_ENZYME_1}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta
			
			mkdir -p ${bionanoPath}/${SC_BIONANO_ENZYME_2}
			
			#---------------------------  enzyme-2
			fname="${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/out/${SC_BIONANO_ENZYME_2}/hybrid_scaffolds_M1/${PROJECT_ID_CAPS}_REFINEFINAL1_bppAdjust_cmap_${REF_NAME}_NGScontigs_HYBRID_SCAFFOLD_NCBI.fasta"
			agp="${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/out/${SC_BIONANO_ENZYME_2}/hybrid_scaffolds_M1/${PROJECT_ID_CAPS}_REFINEFINAL1_bppAdjust_cmap_${REF_NAME}_NGScontigs_HYBRID_SCAFFOLD.agp"
			if [[ ! -f ${fname} ]]
			then
				(>&2 echo "ERROR - Cannot find file ${fname}")
	  			exit 1
			fi
			# enzyme-2: bionano scaffolds
			fext="n"
			cp ${fname} ${bionanoPath}/${SC_BIONANO_ENZYME_2}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta
			gzip -c ${bionanoPath}/${SC_BIONANO_ENZYME_2}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta > ${bionanoPath}/${SC_BIONANO_ENZYME_2}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fa.gz
			cat ${bionanoPath}/${SC_BIONANO_ENZYME_2}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${bionanoPath}/${SC_BIONANO_ENZYME_2}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.stats
			${QUAST_PATH}/quast.py -t 1 -s -e --fast --est-ref-size ${gsize} -o ${bionanoPath}/${SC_BIONANO_ENZYME_2}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset} ${bionanoPath}/${SC_BIONANO_ENZYME_2}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta
			cp ${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/out/${SC_BIONANO_ENZYME_2}/hybrid_scaffolds_M1/hybrid_scaffold_informatics_report.txt ${bionanoPath}/${SC_BIONANO_ENZYME_2}/
			cp ${agp} ${bionanoPath}/${SC_BIONANO_ENZYME_2}
			# enzyme-2: cut bionano conflicts in input assembly
			fext="C" 
			cp ${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/out/${SC_BIONANO_ENZYME_2}/hybrid_scaffolds_M1/$(basename ${SC_BIONANO_REF}).cut.fasta ${bionanoPath}/${SC_BIONANO_ENZYME_2}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta
			cat ${bionanoPath}/${SC_BIONANO_ENZYME_2}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${bionanoPath}/${SC_BIONANO_ENZYME_2}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.stats
			cp ${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/out/${SC_BIONANO_ENZYME_2}/hybrid_scaffolds_M1/conflicts.txt ${bionanoPath}/${SC_BIONANO_ENZYME_2}/
			${QUAST_PATH}/quast.py -t 1 -e --fast --est-ref-size ${gsize} -o ${bionanoPath}/${SC_BIONANO_ENZYME_2}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset} ${bionanoPath}/${SC_BIONANO_ENZYME_2}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta
			# enzyme-2: bionano not scaffolded 
			fname="${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/out/${SC_BIONANO_ENZYME_2}/hybrid_scaffolds_M1/${PROJECT_ID_CAPS}_REFINEFINAL1_bppAdjust_cmap_${REF_NAME}_NGScontigs_HYBRID_SCAFFOLD_NOT_SCAFFOLDED.fasta"
			if [[ ! -f ${fname} ]]
			then
				(>&2 echo "ERROR - Cannot find file ${fname}")
	  			exit 1
			fi
			fext="u" ### bionano not scaffolded
			cp ${fname} ${bionanoPath}/${SC_BIONANO_ENZYME_2}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta
			gzip -c ${bionanoPath}/${SC_BIONANO_ENZYME_2}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta > ${bionanoPath}/${SC_BIONANO_ENZYME_2}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fa.gz
			cat ${bionanoPath}/${SC_BIONANO_ENZYME_2}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${bionanoPath}/${SC_BIONANO_ENZYME_2}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.stats
			${QUAST_PATH}/quast.py -t 1 -s -e --fast --est-ref-size ${gsize} -o ${bionanoPath}/${SC_BIONANO_ENZYME_2}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset} ${bionanoPath}/${SC_BIONANO_ENZYME_2}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta
			# enzyme-2: bionano combined (scaffolded + not scaffolded)
			fext="N" 
			cat ${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/out/${SC_BIONANO_ENZYME_2}/hybrid_scaffolds_M1/${PROJECT_ID_CAPS}_REFINEFINAL1_bppAdjust_cmap_${REF_NAME}_NGScontigs_HYBRID_SCAFFOLD_NCBI.fasta > ${bionanoPath}/${SC_BIONANO_ENZYME_2}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta
			cat ${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/out/${SC_BIONANO_ENZYME_2}/hybrid_scaffolds_M1/${PROJECT_ID_CAPS}_REFINEFINAL1_bppAdjust_cmap_${REF_NAME}_NGScontigs_HYBRID_SCAFFOLD_NOT_SCAFFOLDED.fasta >> ${bionanoPath}/${SC_BIONANO_ENZYME_2}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta
			gzip -c ${bionanoPath}/${SC_BIONANO_ENZYME_2}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta > ${bionanoPath}/${SC_BIONANO_ENZYME_2}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fa.gz
			cat ${bionanoPath}/${SC_BIONANO_ENZYME_2}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${bionanoPath}/${SC_BIONANO_ENZYME_2}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.stats
			${QUAST_PATH}/quast.py -t 1 -s -e --fast --est-ref-size ${gsize} -o ${bionanoPath}/${SC_BIONANO_ENZYME_2}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset} ${bionanoPath}/${SC_BIONANO_ENZYME_2}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta
			
			#--------------------------- both enzymes 
			
			fname="${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/out/two_enzyme_hybrid_scaffold_M1/AGPExport/$(basename ${SC_BIONANO_REF}).cut_${SC_BIONANO_ENZYME_1}_${SC_BIONANO_ENZYME_2}_0kb_0labels_NGS_contigs_HYBRID_Export_NCBI.fasta"
			agp="${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/out/two_enzyme_hybrid_scaffold_M1/AGPExport/$(basename ${SC_BIONANO_REF}).cut_${SC_BIONANO_ENZYME_1}_${SC_BIONANO_ENZYME_2}_0kb_0labels_NGS_contigs_HYBRID_Export.agp"
			if [[ ! -f ${fname} ]]
			then
				(>&2 echo "ERROR - Cannot find file ${fname}")
	  			exit 1
			fi
			# both enzymes: bionano scaffolds
			fext="n"
			cp ${fname} ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta
			gzip -c ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta > ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fa.gz
			cat ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.stats
			${QUAST_PATH}/quast.py -t 1 -s -e --fast --est-ref-size ${gsize} -o ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset} ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta
			cp ${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/out/two_enzyme_hybrid_scaffold_M1/hybrid_scaffold_informatics_report.txt ${bionanoPath}/
			cp ${agp} ${bionanoPath}/
			# both enzymes: cut bionano conflicts in input assembly
			fext="C" 
			cp ${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/out/two_enzyme_hybrid_scaffold_M1/AGPExport/$(basename ${SC_BIONANO_REF}).cut.fasta ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta
			cat ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.stats
			cp ${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/out/two_enzyme_hybrid_scaffold_M1/CombinedConflictsCut/Combined_conflicts_cut_status_${SC_BIONANO_ENZYME_1}.txt ${bionanoPath}/
			cp ${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/out/two_enzyme_hybrid_scaffold_M1/CombinedConflictsCut/Combined_conflicts_cut_status_${SC_BIONANO_ENZYME_2}.txt ${bionanoPath}/
			${QUAST_PATH}/quast.py -t 1 -e --fast --est-ref-size ${gsize} -o ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset} ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta
			# both enzymes: bionano not scaffolded 
			fname="${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/out/two_enzyme_hybrid_scaffold_M1/AGPExport/$(basename ${SC_BIONANO_REF}).cut_${SC_BIONANO_ENZYME_1}_${SC_BIONANO_ENZYME_2}_0kb_0labels_NGS_contigs_HYBRID_Export_NOT_SCAFFOLDED.fasta"
			if [[ ! -f ${fname} ]]
			then
				(>&2 echo "ERROR - Cannot find file ${fname}")
	  			exit 1
			fi
			fext="u" ### bionano not scaffolded
			cp ${fname} ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta
			gzip -c ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta > ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fa.gz
			cat ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.stats
			${QUAST_PATH}/quast.py -t 1 -s -e --fast --est-ref-size ${gsize} -o ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset} ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta
			# both enzymes: bionano combined (scaffolded + not scaffolded)
			fext="N" 
			cat ${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/out/two_enzyme_hybrid_scaffold_M1/AGPExport/$(basename ${SC_BIONANO_REF}).cut_${SC_BIONANO_ENZYME_1}_${SC_BIONANO_ENZYME_2}_0kb_0labels_NGS_contigs_HYBRID_Export_NCBI.fasta > ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta
			cat ${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/out/two_enzyme_hybrid_scaffold_M1/AGPExport/$(basename ${SC_BIONANO_REF}).cut_${SC_BIONANO_ENZYME_1}_${SC_BIONANO_ENZYME_2}_0kb_0labels_NGS_contigs_HYBRID_Export_NOT_SCAFFOLDED.fasta >> ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta
			gzip -c ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta > ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fa.gz
			cat ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.stats
			${QUAST_PATH}/quast.py -t 1 -s -e --fast --est-ref-size ${gsize} -o ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset} ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta
			
		else 
			### single enzyme workflow
			fname="${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/out/hybrid_scaffolds/${PROJECT_ID_CAPS}_REFINEFINAL1_bppAdjust_cmap_${REF_NAME}_NGScontigs_HYBRID_SCAFFOLD_NCBI.fasta"
			if [[ ! -f ${fname} ]]
			then
				(>&2 echo "ERROR - Cannot find file ${fname}")
	  			exit 1
			fi
			
			fext="n"
			# bionano scaffolds
			cat ${fname} > ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta
			gzip -c ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta > ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fa.gz
			cat ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.stats
			${QUAST_PATH}/quast.py -t 1 -s -e --fast --est-ref-size ${gsize} -o ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset} ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta
			cp ${config} ${bionanoPath}/$(date '+%Y-%m-%d_%H-%M-%S')_$(basename ${config})
			cp ${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/out/hybrid_scaffolds/hybrid_scaffold_informatics_report.txt ${bionanoPath}
			cp ${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/out/hybrid_scaffolds/${PROJECT_ID_CAPS}_REFINEFINAL1_bppAdjust_cmap_${REF_NAME}_NGScontigs_HYBRID_SCAFFOLD.agp ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.agp
			
			# cut bionano conflicts in input assembly
			fext="C" 
			cp ${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/out/hybrid_scaffolds/$(basename ${SC_BIONANO_REF}).cut.fasta ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta
			cat ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.stats
			cp ${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/out/hybrid_scaffolds/conflicts.txt ${bionanoPath}
			${QUAST_PATH}/quast.py -t 1 -e --fast --est-ref-size ${gsize} -o ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset} ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta
			
			# bionano not scaffolded 
			fname="${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/out/hybrid_scaffolds/${PROJECT_ID_CAPS}_REFINEFINAL1_bppAdjust_cmap_${REF_NAME}_NGScontigs_HYBRID_SCAFFOLD_NOT_SCAFFOLDED.fasta"
			if [[ ! -f ${fname} ]]
			then
				(>&2 echo "ERROR - Cannot find file ${fname}")
	  			exit 1
			fi
			fext="u" ### bionano not scaffolded
			cat ${fname} > ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta
			gzip -c ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta > ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fa.gz
			cat ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.stats
			${QUAST_PATH}/quast.py -t 1 -s -e --fast --est-ref-size ${gsize} -o ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset} ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta
					
			# bionano combined (scaffolded + not scaffolded)
			fext="N" 
			cat ${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/out/hybrid_scaffolds/${PROJECT_ID_CAPS}_REFINEFINAL1_bppAdjust_cmap_${REF_NAME}_NGScontigs_HYBRID_SCAFFOLD_NCBI.fasta > ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta
			cat ${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/out/hybrid_scaffolds/${PROJECT_ID_CAPS}_REFINEFINAL1_bppAdjust_cmap_${REF_NAME}_NGScontigs_HYBRID_SCAFFOLD_NOT_SCAFFOLDED.fasta >> ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta
			gzip -c ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta > ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fa.gz
			cat ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.stats
			${QUAST_PATH}/quast.py -t 1 -s -e --fast --est-ref-size ${gsize} -o ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset} ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${prevExt}${fext}.${cset}.fasta
		fi
	else
		(>&2 echo "ERROR - directory ${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID} not available")
  		exit 1
	fi	
elif [[ ${phase} -eq 14 ]] ## HIC SALSA scaffolding 
then
	if [[ -d ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID} ]]
	then
		hicSalsaPath=stats/contigs/${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}
		
		prevExt=$(basename ${SC_HIC_REF%.fasta} | awk -F '[_.]' '{print $(NF-1)}')
		cset=$(basename ${SC_HIC_REF%.fasta} | awk -F '[_.]' '{print $(NF)}')
		fext="S"
				
		if [[ -d ${hicSalsaPath} ]]
		then
			mv ${hicSalsaPath} ${hicSalsaPath}_$(date '+%Y-%m-%d_%H-%M-%S')	
		fi
		mkdir -p ${hicSalsaPath}
		
		fname="${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/out/scaffolds_FINAL.fasta"
		if [[ ! -f ${fname} ]]
		then
			(>&2 echo "ERROR - Cannot find file ${fname}")
  			exit 1
		fi
		
		# HIC SALSA scaffolds
		cat ${fname} > ${hicSalsaPath}/${PROJECT_ID}_${SC_HIC_OUTDIR}_${prevExt}${fext}.${cset}.fasta
		gzip -c ${hicSalsaPath}/${PROJECT_ID}_${SC_HIC_OUTDIR}_${prevExt}${fext}.${cset}.fasta > ${hicSalsaPath}/${PROJECT_ID}_${SC_HIC_OUTDIR}_${prevExt}${fext}.${cset}.fa.gz
		cat ${hicSalsaPath}/${PROJECT_ID}_${SC_HIC_OUTDIR}_${prevExt}${fext}.${cset}.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${hicSalsaPath}/${PROJECT_ID}_${SC_HIC_OUTDIR}_${prevExt}${fext}.${cset}.stats
		${QUAST_PATH}/quast.py -o ${hicSalsaPath} -t 1 -s -e --est-ref-size ${gsize} -o ${hicSalsaPath}/${PROJECT_ID}_${SC_HIC_OUTDIR}_${prevExt}${fext}.${cset} ${hicSalsaPath}/${PROJECT_ID}_${SC_HIC_OUTDIR}_${prevExt}${fext}.${cset}.fasta
		cp ${config} ${hicSalsaPath}/$(date '+%Y-%m-%d_%H-%M-%S')_$(basename ${config})
		cp ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/out/input_breaks ${hicSalsaPath}				
	else
		(>&2 echo "ERROR - directory ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID} not available")
  		exit 1
	fi		
elif [[ ${phase} -eq 15 ]] ## PHASING 
then
	if [[ -d ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID} ]]
	then
		
		phasePath=stats/contigs/${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}
		prevExt=$(basename ${CT_PHASE_REFFASTA%.fasta} | awk -F '[_.]' '{print $(NF-1)}')
		cset=$(basename ${CT_PHASE_REFFASTA%.fasta} | awk -F '[_.]' '{print $(NF)}')
		
		if [[ ${CT_PHASE_TYPE} -eq 1 ]]
		then
			fext="s"
			
			if [[ -d ${phasePath} ]]
			then
				mv ${phasePath} ${phasePath}_$(date '+%Y-%m-%d_%H-%M-%S')	
			fi
			mkdir -p ${phasePath}
			
			fname1="${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/${PROJECT_ID}_phased.h1.fasta"
			fname2="${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/${PROJECT_ID}_phased.h2.fasta"
			if [[ ! -f ${fname1} || ! -f ${fname2} ]]
			then
				(>&2 echo "ERROR - Cannot find phased files ${fname1} and / or ${fname2}")
	  			exit 1
			fi
			
			cp ${fname1} ${phasePath}/${PROJECT_ID}_${CT_PHASE_OUTDIR}_${prevExt}${fext}.${cset}1.fasta
			gzip -c ${phasePath}/${PROJECT_ID}_${CT_PHASE_OUTDIR}_${prevExt}${fext}.${cset}1.fasta > ${phasePath}/${PROJECT_ID}_${CT_PHASE_OUTDIR}_${prevExt}${fext}.${cset}1.fa.gz	
			cat ${phasePath}/${PROJECT_ID}_${CT_PHASE_OUTDIR}_${prevExt}${fext}.${cset}1.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${phasePath}/${PROJECT_ID}_${CT_PHASE_OUTDIR}_${prevExt}${fext}.${cset}1.stats
			${QUAST_PATH}/quast.py -t 1 -s -e --est-ref-size ${gsize} -o ${phasePath}/${PROJECT_ID}_${CT_PHASE_OUTDIR}_${prevExt}${fext}.${cset}1 ${phasePath}/${PROJECT_ID}_${CT_PHASE_OUTDIR}_${prevExt}${fext}.${cset}1.fasta

			cp ${fname2} ${phasePath}/${PROJECT_ID}_${CT_PHASE_OUTDIR}_${prevExt}${fext}.${cset}2.fasta
			gzip -c ${phasePath}/${PROJECT_ID}_${CT_PHASE_OUTDIR}_${prevExt}${fext}.${cset}2.fasta > ${phasePath}/${PROJECT_ID}_${CT_PHASE_OUTDIR}_${prevExt}${fext}.${cset}2.fa.gz	
			cat ${phasePath}/${PROJECT_ID}_${CT_PHASE_OUTDIR}_${prevExt}${fext}.${cset}2.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${phasePath}/${PROJECT_ID}_${CT_PHASE_OUTDIR}_${prevExt}${fext}.${cset}2.stats
			${QUAST_PATH}/quast.py -t 1 -s -e --est-ref-size ${gsize} -o ${phasePath}/${PROJECT_ID}_${CT_PHASE_OUTDIR}_${prevExt}${fext}.${cset}2 ${phasePath}/${PROJECT_ID}_${CT_PHASE_OUTDIR}_${prevExt}${fext}.${cset}2.fasta
		
			cp ${config} ${phasePath}/$(date '+%Y-%m-%d_%H-%M-%S')_$(basename ${config})
			
			## add phasing stats
			awk -F , '{ for (i=1; i<=NF; i++)  { a[NR,i] = $i } } NF>p { p = NF } END { for(j=1; j<=p; j++) { str=a[1,j]; for(i=2; i<=NR; i++){ str=str"\t"a[i,j]; } print str } }' ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/bams/10x_${PROJECT_ID}_longrangerWgs/outs/summary.csv >  ${phasePath}/${PROJECT_ID}_${CT_PHASE_OUTDIR}_${prevExt}${fext}.${cset}.summary
			## check if some contigs were ignored				
			if [[ -s ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/ref/ignore.fasta ]]
			then
				echo "Unphased contigs num: $(grep -c -e \">\" ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/ref/ignore.fasta) size: $(du -h ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID}/ref/ignore.fasta | awk '{print $1}')" > ${phasePath}/${PROJECT_ID}_${CT_PHASE_OUTDIR}_${prevExt}${fext}.${cset}.error
			fi			
		else
			(>&2 echo "ERROR - CT_PHASE_TYPE: ${CT_PHASE_TYPE} not supported yet!")
  			exit 1
		fi
	else
		(>&2 echo "ERROR - directory ${CT_PHASE_OUTDIR}/phase_${CT_PHASE_RUNID} not available")
  		exit 1
	fi	
else
	(>&2 echo "Unknow Phase: ${phase}")
	exit 1	
fi