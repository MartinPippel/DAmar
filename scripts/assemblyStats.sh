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
			
			cat ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${name}${pathID}${arrowExtension}/ALL_${name}${pathID}${arrowExtension}.arrow.fa	| sed -e "s:[|_]arrow::g"					        		
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
			
			cat ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${name}${pathID}${arrowExtension}/ALL_${name}${pathID}${arrowExtension}.arrow.fa	| sed -e "s:[|_]arrow::g"					        		
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
			
			cat ${PB_ARROW_OUTDIR}/arrow_${PB_ARROW_RUNID}/${name}${pathID}${arrowExtension}/ALL_${name}${pathID}${arrowExtension}.arrow.fa	| sed -e "s:[|_]arrow::g"					        		
		done > ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.fasta	
        gzip -c ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.fasta > ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.fa.gz
        if [[ -s ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.fasta ]]
        then 
        	cat ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.stats
        	allBases=$(cat ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.fasta | awk '{ if ($1 ~ /^>/) ; else print $0;}' | tr -d "\n" | wc -m)
        	arrowedBases=$(cat ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.fasta | awk '{ if ($1 ~ /^>/) ; else print $0;}' | tr -d "\nacgt" | wc -m)
        	frac=$(echo "scale=2;${arrowedBases}*100/\${allBases}" | bc)
        	echo "allBases ${allBases} arrowedBases ${arrowedBases} arrowedFraction ${frac}%" >> ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.stats
        	${SUBMIT_SCRIPTS_PATH}/trimLowerCaseTips.py ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.fasta ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.T.fasta ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.T.failed.fasta
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
        	${SUBMIT_SCRIPTS_PATH}/trimLowerCaseTips.py ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.p.fasta ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.pT.fasta ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.pT.failed.fasta
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
        	${SUBMIT_SCRIPTS_PATH}/trimLowerCaseTips.py ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.a.fasta ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.aT.fasta ${arrowPath}/${PROJECT_ID}_${FIX_FILT_OUTDIR}_${fext}.aT.failed.fasta
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
	if [[ -d ${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID} ]]
	then
		freebayesPath=stats/contigs/${CT_FREEBAYES_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID}
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
		
		cat ${fname} > ${freebayesPath}/${PROJECT_ID}_${CT_FREEBAYES_OUTDIR}_${fext}.p.fasta
		gzip -c ${freebayesPath}/${PROJECT_ID}_${CT_FREEBAYES_OUTDIR}_${fext}.p.fasta > ${freebayesPath}/${PROJECT_ID}_${CT_FREEBAYES_OUTDIR}_${fext}.p.fa.gz
		cat ${freebayesPath}/${PROJECT_ID}_${CT_FREEBAYES_OUTDIR}_${fext}.p.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${freebayesPath}/${PROJECT_ID}_${CT_FREEBAYES_OUTDIR}_${fext}.p.stats
		cp ${config} ${freebayesPath}/$(date '+%Y-%m-%d_%H-%M-%S')_$(basename ${config})
	else
		(>&2 echo "ERROR - directory ${FIX_FILT_OUTDIR}/freebayes_${CT_FREEBAYES_RUNID} not available")
  		exit 1
	fi
elif [[ ${phase} -eq 12 ]] ## scaff10x scaffolding 
then
	if [[ -d ${SC_SCAFF10X_OUTDIR}/scaff10x_${SC_SCAFF10X_RUNID} ]]
	then
		scaff10xPath=stats/contigs/${SC_SCAFF10X_OUTDIR}/scaff10x_${SC_SCAFF10X_RUNID}
		fext="x"
		
		mkdir -p ${scaff10xPath}
		
		fname="${SC_SCAFF10X_OUTDIR}/scaff10x_${SC_SCAFF10X_RUNID}/${PROJECT_ID}_${SC_SCAFF10X_OUTDIR}_${fext}.p.fasta"
		if [[ ! -f ${fname} ]]
		then
			(>&2 echo "ERROR - Cannot find file ${fname}")
  			exit 1
		fi
		
		# scaff10x scaffolds
		cat ${fname} > ${scaff10xPath}/${PROJECT_ID}_${SC_SCAFF10X_OUTDIR}_${fext}.p.fasta
		gzip -c ${scaff10xPath}/${PROJECT_ID}_${SC_SCAFF10X_OUTDIR}_${fext}.p.fasta > ${scaff10xPath}/${PROJECT_ID}_${SC_SCAFF10X_OUTDIR}_${fext}.p.fa.gz
		cat ${scaff10xPath}/${PROJECT_ID}_${SC_SCAFF10X_OUTDIR}_${fext}.p.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${scaff10xPath}/${PROJECT_ID}_${SC_SCAFF10X_OUTDIR}_${fext}.p.stats
		${QUAST_PATH}/quast.py -o ${scaff10xPath} -t 1 -s -e --est-ref-size ${gsize} ${scaff10xPath}/${PROJECT_ID}_${SC_SCAFF10X_OUTDIR}_${fext}.p.fasta
		cp ${config} ${scaff10xPath}/$(date '+%Y-%m-%d_%H-%M-%S')_$(basename ${config})
		# break10x based in scaff10x input
		if [[ -f ${SC_SCAFF10X_OUTDIR}/scaff10x_${SC_SCAFF10X_RUNID}/$(basename ${SC_SCAFF10X_REF%.*})_break10x.fasta && -f ${SC_SCAFF10X_OUTDIR}/scaff10x_${SC_SCAFF10X_RUNID}/$(basename ${SC_SCAFF10X_REF%.*})_break10x.breaks ]]
		then
			cp ${SC_SCAFF10X_OUTDIR}/scaff10x_${SC_SCAFF10X_RUNID}/$(basename ${SC_SCAFF10X_REF%.*})_break10x.fasta ${scaff10xPath}/${PROJECT_ID}_${SC_SCAFF10X_OUTDIR}_b${fext}.p.fasta
			cp ${SC_SCAFF10X_OUTDIR}/scaff10x_${SC_SCAFF10X_RUNID}/$(basename ${SC_SCAFF10X_REF%.*})_break10x.breaks ${scaff10xPath}/${PROJECT_ID}_${SC_SCAFF10X_OUTDIR}_b${fext}.p.breaks
			gzip -c ${scaff10xPath}/${PROJECT_ID}_${SC_SCAFF10X_OUTDIR}_b${fext}.p.fasta > ${scaff10xPath}/${PROJECT_ID}_${SC_SCAFF10X_OUTDIR}_b${fext}.p.fa.gz
			cat ${scaff10xPath}/${PROJECT_ID}_${SC_SCAFF10X_OUTDIR}_b${fext}.p.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${scaff10xPath}/${PROJECT_ID}_${SC_SCAFF10X_OUTDIR}_b${fext}.p.stats
			${QUAST_PATH}/quast.py -o ${scaff10xPath} -t 1 -s -e --est-ref-size ${gsize} ${scaff10xPath}/${PROJECT_ID}_${SC_SCAFF10X_OUTDIR}_b${fext}.p.fasta
		fi
		# break10x based on scaff10x output
		if [[ -f ${SC_SCAFF10X_OUTDIR}/scaff10x_${SC_SCAFF10X_RUNID}/${PROJECT_ID}_${SC_SCAFF10X_OUTDIR}_x.p_break10x.fasta && -f ${SC_SCAFF10X_OUTDIR}/scaff10x_${SC_SCAFF10X_RUNID}/${PROJECT_ID}_${SC_SCAFF10X_OUTDIR}_x.p_break10x.breaks ]]
		then
			cp ${SC_SCAFF10X_OUTDIR}/scaff10x_${SC_SCAFF10X_RUNID}/${PROJECT_ID}_${SC_SCAFF10X_OUTDIR}_x.p_break10x.fasta ${scaff10xPath}/${PROJECT_ID}_${SC_SCAFF10X_OUTDIR}_B${fext}.p.fasta
			cp ${SC_SCAFF10X_OUTDIR}/scaff10x_${SC_SCAFF10X_RUNID}/${PROJECT_ID}_${SC_SCAFF10X_OUTDIR}_x.p_break10x.breaks ${scaff10xPath}/${PROJECT_ID}_${SC_SCAFF10X_OUTDIR}_B${fext}.p.breaks
			gzip -c ${scaff10xPath}/${PROJECT_ID}_${SC_SCAFF10X_OUTDIR}_B${fext}.p.fasta > ${scaff10xPath}/${PROJECT_ID}_${SC_SCAFF10X_OUTDIR}_B${fext}.p.fa.gz
			cat ${scaff10xPath}/${PROJECT_ID}_${SC_SCAFF10X_OUTDIR}_B${fext}.p.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${scaff10xPath}/${PROJECT_ID}_${SC_SCAFF10X_OUTDIR}_B${fext}.p.stats
			${QUAST_PATH}/quast.py -o ${scaff10xPath} -t 1 -s -e --est-ref-size ${gsize} ${scaff10xPath}/${PROJECT_ID}_${SC_SCAFF10X_OUTDIR}_B${fext}.p.fasta				
		fi 		
	else
		(>&2 echo "ERROR - directory ${SC_SCAFF10X_OUTDIR}/scaff10x_${SC_SCAFF10X_RUNID} not available")
  		exit 1
	fi	
elif [[ ${phase} -eq 13 ]] ## bionano scaffolding 
then
	if [[ -d ${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID} ]]
	then
		bionanoPath=stats/contigs/${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}
		
		fext="b"
				
		mkdir -p ${bionanoPath}
		
		PROJECT_ID_CAPS=$(echo ${PROJECT_ID} | awk '{print toupper($0)}')
		REF_NAME=$(basename ${SC_BIONANO_REF} | tr '.' '_')
		fname="${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/out/hybrid_scaffolds/${PROJECT_ID_CAPS}_REFINEFINAL1_bppAdjust_cmap_${REF_NAME}_NGScontigs_HYBRID_SCAFFOLD.fasta"
		if [[ ! -f ${fname} ]]
		then
			(>&2 echo "ERROR - Cannot find file ${fname}")
  			exit 1
		fi
		
		# bionano scaffolds
		cat ${fname} > ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_s${fext}.p.fasta
		gzip -c ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_s${fext}.p.fasta > ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_s${fext}.p.fa.gz
		cat ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_s${fext}.p.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_s${fext}.p.stats
		${QUAST_PATH}/quast.py -o ${bionanoPath} -t 1 -s -e --est-ref-size ${gsize} ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_s${fext}.p.fasta
		cp ${config} ${bionanoPath}/$(date '+%Y-%m-%d_%H-%M-%S')_$(basename ${config})
		cp ${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/out/hybrid_scaffolds/hybrid_scaffold_informatics_report.txt ${bionanoPath}
		cp ${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/out/hybrid_scaffolds/hybrid_scaffold_informatics_report.txt ${bionanoPath}
		cp ${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/out/hybrid_scaffolds/$(basename ${SC_BIONANO_REF}).cut.fasta ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_c${fext}.p.fasta
		cp ${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/out/hybrid_scaffolds/conflicts.txt ${bionanoPath}
		
		# bionano not scaffolded 
		fname="${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/out/hybrid_scaffolds/${PROJECT_ID_CAPS}_REFINEFINAL1_bppAdjust_cmap_${REF_NAME}_NGScontigs_HYBRID_SCAFFOLD_NOT_SCAFFOLDED.fasta"
		if [[ ! -f ${fname} ]]
		then
			(>&2 echo "ERROR - Cannot find file ${fname}")
  			exit 1
		fi
		cat ${fname} > ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_n${fext}.p.fasta
		gzip -c ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_n${fext}.p.fasta > ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_n${fext}.p.fa.gz
		cat ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_n${fext}.p.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_n${fext}.p.stats
		${QUAST_PATH}/quast.py -o ${bionanoPath} -t 1 -s -e --est-ref-size ${gsize} ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_n${fext}.p.fasta
				
		# bionano combined (scaffolded + not scaffolded)
		cat ${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/out/hybrid_scaffolds/${PROJECT_ID_CAPS}_REFINEFINAL1_bppAdjust_cmap_${REF_NAME}_NGScontigs_HYBRID_SCAFFOLD.fasta > ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${fext}.p.fasta
		cat ${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID}/out/hybrid_scaffolds/${PROJECT_ID_CAPS}_REFINEFINAL1_bppAdjust_cmap_${REF_NAME}_NGScontigs_HYBRID_SCAFFOLD_NOT_SCAFFOLDED.fasta >> ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${fext}.p.fasta
		gzip -c ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${fext}.p.fasta > ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${fext}.p.fa.gz
		cat ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${fext}.p.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${fext}.p.stats
		${QUAST_PATH}/quast.py -o ${bionanoPath} -t 1 -s -e --est-ref-size ${gsize} ${bionanoPath}/${PROJECT_ID}_${SC_BIONANO_OUTDIR}_${fext}.p.fasta
	else
		(>&2 echo "ERROR - directory ${SC_BIONANO_OUTDIR}/bionano_${SC_BIONANO_RUNID} not available")
  		exit 1
	fi	
elif [[ ${phase} -eq 14 ]] ## HIC SALSA scaffolding 
then
	if [[ -d ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID} ]]
	then
		hicSalsaPath=stats/contigs/${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}
		
		fext="sH"
				
		mkdir -p ${hicSalsaPath}
		
		REF_NAME=$(basename ${SC_HIC_REFFASTA})
		fname="${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/${REF_NAME}/scaffolds_FINAL.fasta"
		if [[ ! -f ${fname} ]]
		then
			(>&2 echo "ERROR - Cannot find file ${fname}")
  			exit 1
		fi
		
		# HIC SALSA scaffolds
		cat ${fname} > ${hicSalsaPath}/${PROJECT_ID}_${SC_HIC_OUTDIR}_${fext}.p.fasta
		gzip -c ${hicSalsaPath}/${PROJECT_ID}_${SC_HIC_OUTDIR}_${fext}.p.fasta > ${hicSalsaPath}/${PROJECT_ID}_${SC_HIC_OUTDIR}_${fext}.p.fa.gz
		cat ${hicSalsaPath}/${PROJECT_ID}_${SC_HIC_OUTDIR}_${fext}.p.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > ${hicSalsaPath}/${PROJECT_ID}_${SC_HIC_OUTDIR}_${fext}.p.stats
		${QUAST_PATH}/quast.py -o ${hicSalsaPath} -t 1 -s -e --est-ref-size ${gsize} ${hicSalsaPath}/${PROJECT_ID}_${SC_HIC_OUTDIR}_${fext}.p.fasta
		cp ${config} ${hicSalsaPath}/$(date '+%Y-%m-%d_%H-%M-%S')_$(basename ${config})
		cp ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID}/${REF_NAME}/input_breaks ${hicSalsaPath}				
	else
		(>&2 echo "ERROR - directory ${SC_HIC_OUTDIR}/hic_${SC_HIC_RUNID} not available")
  		exit 1
	fi			
else
	(>&2 echo "Unknow Phase: ${phase}")
	exit 1	
fi
			
	
