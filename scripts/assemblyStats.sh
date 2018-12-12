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
	if [[ -d ${FIX_FILT_OUTDIR}/correction/contigs ]]
	then
		p=stats/contigs/${FIX_FILT_OUTDIR}/corr
		mkdir -p ${p}
		
		for y in ${FIX_FILT_OUTDIR}/correction/contigs/*.fasta
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
		(>&2 echo "ERROR - directory ${FIX_FILT_OUTDIR}/correction/contigs not available")
  		exit 1
	fi
	
elif [[ ${phase} -eq 8 ]] ## marvel corrected assembly stats after contig analysis  (last step after CTanalysis)
then
	echo "#todo" 
elif [[ ${phase} -eq 9 ]] ## assembly stats after PacBio Arrow Correction 
then
	if [[ -d ${FIX_FILT_OUTDIR}/arrow_${PB_ARROW_RUNID} ]]
	then
		
		rawPath=stats/contigs/${FIX_FILT_OUTDIR}/raw
		arrowPath=stats/contigs/${FIX_FILT_OUTDIR}/arrow_${PB_ARROW_RUNID}
		fext="a"
		if [[ ${PB_ARROW_RUNID} -eq 2 ]]
		then 
			fext="A"
		elif [[ ${PB_ARROW_RUNID} -gt 2 ]]
		then 
			fext="@"		
		fi
		
		if [[ ! -d ${rawPath} || ! -f ${rawPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_r.p.fasta || ! -f ${rawPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_r.a.fasta ]]
        then 
	     	(>&2 echo "ERROR - stats folder or assembly staticstics are missing. Run last step of touring first.")
    	    exit 1        			
        fi
        
        mkdir -p ${arrowPath}        

		## primary 
		for z in $(cat ${rawPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_r.p.header)
		do
			name=$(echo ${z} | awk -F \_ '{$NF=""; OFS="_"; print $0}')
			pathID=$(echo ${z} | awk -F \_ '{print $NF}')
			
			arrowExtension=""
			a=1
			while [[ $a -lt ${PB_ARROW_RUNID} ]]
			do
				arrowExtension="${arrowExtension}_arrow"
				a=$(($a+1))	
			done
			
			if [[ ! -d ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${name}${pathID}${arrowExtension} ]]
			then 
				echo "WARNING - Missing directory ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${name}${pathID}${arrowExtension}." >> ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.p.elog
				continue
			fi  
			
			if [[ ! -f ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${name}${pathID}${arrowExtension}/ALL_${name}${pathID}${arrowExtension}.arrow.fa ]]
			then 
			 	echo "WARNING - Arrow failed. Missing file ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${name}${pathID}${arrowExtension}/ALL_${name}${pathID}${arrowExtension}.arrow.fa." >> ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.p.elog
				continue
			fi  
			
			cat ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${name}${pathID}${arrowExtension}/ALL_${name}${pathID}${arrowExtension}.arrow.fa	| sed -e "s:|arrow::g"					        		
		done > ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.p.fasta
		gzip -c ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.p.fasta > ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.p.fa.gz
		## alternate
		for z in $(cat ${rawPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_r.a.header)
		do
			name=$(echo ${z} | awk -F \_ '{$NF=""; OFS="_"; print $0}')
			pathID=$(echo ${z} | awk -F \_ '{print $NF}')
			
			arrowExtension=""
			a=1
			while [[ $a -lt ${PB_ARROW_RUNID} ]]
			do
				arrowExtension="${arrowExtension}_arrow"
				a=$(($a+1))	
			done
			
			if [[ ! -d ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${name}${pathID}${arrowExtension} ]]
			then 
				echo "WARNING - Missing directory ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${name}${pathID}${arrowExtension}." >> ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.a.elog
				continue
			fi  
			
			if [[ ! -f ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${name}${pathID}${arrowExtension}/ALL_${name}${pathID}${arrowExtension}.arrow.fa ]]
			then 
			 	echo "WARNING - Arrow failed. Missing file ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${name}${pathID}${arrowExtension}/ALL_${name}${pathID}${arrowExtension}.arrow.fa." >> ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.a.elog
				continue
			fi  
			
			cat ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${name}${pathID}${arrowExtension}/ALL_${name}${pathID}${arrowExtension}.arrow.fa	| sed -e "s:|arrow::g"					        		
		done > ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.a.fasta
		gzip -c ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.a.fasta > ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.a.fa.gz
		## all 	
        for z in $(cat ${rawPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_r.header)
		do
			name=$(echo ${z} | awk -F \_ '{$NF=""; OFS="_"; print $0}')
			pathID=$(echo ${z} | awk -F \_ '{print $NF}')
			
			arrowExtension=""
			a=1
			while [[ $a -lt ${PB_ARROW_RUNID} ]]
			do
				arrowExtension="${arrowExtension}_arrow"
				a=$(($a+1))	
			done
			
			if [[ ! -d ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${name}${pathID}${arrowExtension} ]]
			then 
				echo "WARNING - Missing directory ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${name}${pathID}${arrowExtension}." >> ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.elog
				continue
			fi  
			
			if [[ ! -f ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${name}${pathID}${arrowExtension}/ALL_${name}${pathID}${arrowExtension}.arrow.fa ]]
			then 
			 	echo "WARNING - Arrow failed. Missing file ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${name}${pathID}${arrowExtension}/ALL_${name}${pathID}${arrowExtension}.arrow.fa." >> ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.elog
				continue
			fi  
			
			cat ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${name}${pathID}${arrowExtension}/ALL_${name}${pathID}${arrowExtension}.arrow.fa	| sed -e "s:|arrow::g"					        		
		done > ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.fasta	
        gzip -c ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.fasta > ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.fa.gz
        if [[ -s ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.fasta ]]
        then 
        	cat ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.stats
        	allBases=$(cat ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.fasta | awk '{ if ($1 ~ /^>/) ; else print $0;}' | tr -d "\n" | wc -m)
        	arrowedBases=$(cat ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.fasta | awk '{ if ($1 ~ /^>/) ; else print $0;}' | tr -d "\nacgt" | wc -m)
        	frac=$(echo "scale=2;${arrowedBases}*100/\${allBases}" | bc)
        	echo "allBases ${allBases} arrowedBases ${arrowedBases} arrowedFraction ${frac}%" >> ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.stats
        	${SUBMIT_SCRIPTS_PATH}/trimLowercaseTips.py ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.fasta ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.T.fasta ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.T.failed
        	cat ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.T.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.T.stats
			allBases=$(cat ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.T.fasta | awk '{ if ($1 ~ /^>/) ; else print $0;}' | tr -d "\n" | wc -m)
        	arrowedBases=$(cat ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.T.fasta | awk '{ if ($1 ~ /^>/) ; else print $0;}' | tr -d "\nacgt" | wc -m)
        	frac=$(echo "scale=2;${arrowedBases}*100/${allBases}" | bc)
        	echo "allBases ${allBases} arrowedBases ${arrowedBases} arrowedFraction ${frac}%" >> ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.T.stats
    	fi
        	
        if [[ -s ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.p.fasta ]]
        then
        	cat ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.p.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.p.stats
        	allBases=$(cat ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.p.fasta | awk '{ if ($1 ~ /^>/) ; else print $0;}' | tr -d "\n" | wc -m)
        	arrowedBases=$(cat ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.p.fasta | awk '{ if ($1 ~ /^>/) ; else print $0;}' | tr -d "\nacgt" | wc -m)
        	frac=$(echo "scale=2;${arrowedBases}*100/${allBases}" | bc)
        	echo "allBases ${allBases} arrowedBases ${arrowedBases} arrowedFraction ${frac}%" >> ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.p.stats
        	${SUBMIT_SCRIPTS_PATH}/trimLowercaseTips.py ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.p.fasta ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.pT.fasta ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.pT.failed
        	cat ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.pT.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.pT.stats
			allBases=$(cat ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.pT.fasta | awk '{ if ($1 ~ /^>/) ; else print $0;}' | tr -d "\n" | wc -m)
        	arrowedBases=$(cat ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.pT.fasta | awk '{ if ($1 ~ /^>/) ; else print $0;}' | tr -d "\nacgt" | wc -m)
        	frac=$(echo "scale=2;${arrowedBases}*100/${allBases}" | bc)
        	echo "allBases ${allBases} arrowedBases ${arrowedBases} arrowedFraction ${frac}%" >> ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.pT.stats
    	fi 
       	
        if [[ -s ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.a.fasta ]]
        then
        	cat ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.a.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.a.stats
        	allBases=$(cat ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.a.fasta | awk '{ if ($1 ~ /^>/) ; else print $0;}' | tr -d "\n" | wc -m)
        	arrowedBases=$(cat ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.a.fasta | awk '{ if ($1 ~ /^>/) ; else print $0;}' | tr -d "\nacgt" | wc -m)
        	frac=$(echo "scale=2;${arrowedBases}*100/${allBases}" | bc )
        	echo "allBases ${allBases} arrowedBases ${arrowedBases} arrowedFraction ${frac}%" >> ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.a.stats
        	${SUBMIT_SCRIPTS_PATH}/trimLowercaseTips.py ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.a.fasta ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.aT.fasta ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.aT.failed
        	cat ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.aT.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.aT.stats
			allBases=$(cat ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.aT.fasta | awk '{ if ($1 ~ /^>/) ; else print $0;}' | tr -d "\n" | wc -m)
        	arrowedBases=$(cat ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.aT.fasta | awk '{ if ($1 ~ /^>/) ; else print $0;}' | tr -d "\nacgt" | wc -m)
        	frac=$(echo "scale=2;${arrowedBases}*100/${allBases}" | bc)        	
        	echo "allBases ${allBases} arrowedBases ${arrowedBases} arrowedFraction ${frac}%" >> ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.aT.stats
    	fi	 
	else
		(>&2 echo "ERROR - directory ${FIX_FILT_OUTDIR}/arrow_${PB_ARROW_RUNID} not available")
  		exit 1
	fi
elif [[ ${phase} -eq 11 ]] ## after freebayes polishing 
then
	if [[ -d ${FIX_FILT_OUTDIR}/arrow_${PB_ARROW_RUNID} ]]
	then
		
		rawPath=stats/contigs/${FIX_FILT_OUTDIR}/raw
		freebayesPath=stats/contigs/${FIX_FILT_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}
		fext="f"
		
		mkdir -p ${freebayesPath}
		
		fname=""
		if [[ ${CT_FREEBAYES_TYPE} -eq 0 ]]
		then 
			fname="${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/${PROJECT_ID}_10x.fasta"
		elif [[ ${CT_FREEBAYES_TYPE} -eq 1 ]]
		then 
			fname="${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}/${PROJECT_ID}_hic.fasta"
		fi
		
		if [[ ! -f ${fname} ]]
		then
			(>&2 echo "ERROR - Cannot find polished file ${fname}")
  			exit 1
		fi
		
		cat ${fname} > ${freebayesPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.p.fasta
		gzip -c ${freebayesPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.p.fasta > ${freebayesPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.p.fa.gz
		cat ${freebayesPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.p.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${freebayesPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.p.stats
	else
		(>&2 echo "ERROR - directory ${FIX_FILT_OUTDIR}/arrow_${PB_ARROW_RUNID} not available")
  		exit 1
	fi
			
else
	(>&2 echo "Unknow Phase: ${phase}")
	exit 1	
fi
			
	
