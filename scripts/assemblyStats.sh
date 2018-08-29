#!/bin/bash 

config=$1

if [[ ! -f ${config} ]]
then 
  echo "config ${config} not available"
  exit 1
fi

source ${config}

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

for x in ${FIX_FILT_OUTDIR}_dalign ${FIX_FILT_OUTDIR}_repcomp ${FIX_FILT_OUTDIR}_forcealign
do
	if [[ -d ${x}/tour ]]
	then
		mkdir -p stats/contigs/${x}
		
		cat	${x}/tour/*.fasta > stats/contigs/${x}/${x}.fasta
		${SUBMIT_SCRIPTS_PATH}/splitDiploidAssembly.py ${x} ${gsize} stats/contigs/${x} stats/contigs/${x}/${x}.fasta  
		cat stats/contigs/${x}/${x}.haploid.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > stats/contigs/${x}/${x}.haploid.stats
		cat stats/contigs/${x}/${x}.bubbles.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > stats/contigs/${x}/${x}.bubbles.stats
		cat stats/contigs/${x}/${x}.spurs.fasta | ${SUBMIT_SCRIPTS_PATH}/n50.py ${gsize} > stats/contigs/${x}/${x}.spurs.stats
		
		cp $config stats/contigs/${x}/
	fi	
done 
