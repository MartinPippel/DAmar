#!/bin/bash 

config=$1

if [[ ! -f ${config} ]]
then 
  echo "config ${config} not available"
  exit 1
fi

source ${config}

mkdir -p stats

function createSlurmStats()
{
    slurm_id=$1

    runTable=$(sacct -n -j ${slurm_id} -X -P --format="JobID%30,JobName,Partition,State,NCPUS,ElapsedRaw,CPUTimeRAW,Timelimit,ReqMem,Submit,End")        

        if [[ -z ${runTable} ]]
        then
            (>&2 echo "Unable to get Slurm jobs stats for jobid $slurm_id")
            exit 1
        fi

    ## extract things from runtable
    numJobs=$(echo "${runTable}" | wc -l)
    failJobs=$(echo "${runTable}" | grep -c -v -e COMPLETED)
    CumElapsedRaw=$(echo "${runTable}" | awk -F \| '{s+=$6} END {printf "%d|%.1f|%.1f\n", s, s/60, s/3600}')
    CumCPUTimeRAW=$(echo "${runTable}" | awk -F \| '{s+=$7} END {printf "%d|%.1f|%.1f\n", s, s/60, s/3600}')
	partition=$(echo "${runTable}" |awk -F \| '{if($3 == "gpu") gpu+=1; else if ($3 == "batch") batch+=1; else if ($3 == "long") long+=1; else if ($3 == "bigmem") bigmem+=1; else other+=1} END {printf "%d|%d|%d|%d|%d\n", batch, bigmem, gpu, long, other}')
    ReqMem=$(echo "${runTable}" | head -n1 | awk -F \| '{print $9}')
    NCPUS=$(echo "${runTable}" | head -n1 | awk -F \| '{print $5}')
    Submit=$(echo "${runTable}" | awk -F \| '{print $10}' | xargs -i date -d {} +%s | sort -n | head -n1)
    End=$(echo "${runTable}" | awk -F \| '{print $11}' | xargs -i date -d {} +%s | sort -n -r | head -n1)

    ## extract things from intermediate results - i.e. it only wroking for completed jobs
    runTableIntermediate=$(sacct -n -j ${slurm_id} -p --unit M --format="JobID%30,JobName,AveVMSize,AveRSS,NTasks,AveCPUFreq,ConsumedEnergy,AveDiskRead,AveDiskWrite" | awk -F \| '{if ($2=="batch") print $0}' | sed -e "s:M|:|:g")
    AveVMSize=$(echo "$runTableIntermediate" | awk -F \| '{ if($3 !=0) print $3}'  | awk -F \| -v min=99999999 -v max=0 -v count=0 '{ if ($1 < min) min=$1; if ($1 > max) max=$1; count+=1; avg+=$1;} END {if(count > 0) printf "%4.1f|%4.1f|%4.1f", min, max, avg/count}')
    if [[ -z ${AveVMSize} ]]
    then
        AveVMSize="NotAvailable"
    fi    
    AveRSS=$(echo "$runTableIntermediate" | awk -F \| '{ if($4 !=0) print $4}'  | awk -F \| -v min=99999999 -v max=0 -v count=0 '{ if ($1 < min) min=$1; if ($1 > max) max=$1; count+=1; avg+=$1;} END {if(count > 0) printf "%4.1f|%4.1f|%4.1f", min, max, avg/count}')
    if [[ -z ${AveRSS} ]]
    then
        AveRSS="NotAvailable"
    fi
    NTasks=$(echo "$runTableIntermediate" | awk -F \| '{ if($5 !=0) print $5}'  | awk -F \| -v min=99999999 -v max=0 -v count=0 '{ if ($1 < min) min=$1; if ($1 > max) max=$1; count+=1; avg+=$1;} END {if(count > 0) printf "%4.1f|%4.1f|%4.1f", min, max, avg/count}')
    if [[ -z ${NTasks} ]]
    then
        NTasks="NotAvailable"
    fi
    AveCPUFreq=$(echo "$runTableIntermediate" | awk -F \| '{ if($6 !=0) print $6}'  | awk -F \| -v min=99999999 -v max=0 -v count=0 '{ if ($1 < min) min=$1; if ($1 > max) max=$1; count+=1; avg+=$1;} END {if(count > 0) printf "%4.1f|%4.1f|%4.1f", min/1000, max/1000, avg/count/1000}')
    if [[ -z ${AveCPUFreq} ]]
    then
        AveCPUFreq="NotAvailable"
    fi
    ConsumedEnergy=$(echo "$runTableIntermediate" | awk -F \| '{ if($7 !=0) print $7}'  | awk -F \| -v min=99999999 -v max=0 -v count=0 '{ if ($1 < min) min=$1; if ($1 > max) max=$1; count+=1; avg+=$1;} END {if(count > 0) printf "%4.1f|%4.1f|%4.1f", min, max, avg/count}')
    if [[ -z ${ConsumedEnergy} ]]
    then 
        ConsumedEnergy="NotAvailable"
    fi
    AveDiskRead=$(echo "$runTableIntermediate" | awk -F \| '{ if($8 !=0) print $8}'  | awk -F \| -v min=99999999 -v max=0 -v count=0 '{ if ($1 < min) min=$1; if ($1 > max) max=$1; count+=1; avg+=$1;} END {if(count > 0) printf "%4.1f|%4.1f|%4.1f", min, max, avg/count}')
    if [[ -z ${AveDiskRead} ]]
    then
        AveDiskRead="NotAvailable"
    fi
    AveDiskWrite=$(echo "$runTableIntermediate" | awk -F \| '{ if($9 !=0) print $9}'  | awk -F \| -v min=99999999 -v max=0 -v count=0 '{ if ($1 < min) min=$1; if ($1 > max) max=$1; count+=1; avg+=$1;} END {if(count > 0) printf "%4.1f|%4.1f|%4.1f", min, max, avg/count}')
    if [[ -z ${AveDiskWrite} ]]
    then
        AveDiskWrite="NotAvailable"
    fi

    printf "%15s %25s         %s\n" "SlurmID" 	"${slurm_id}" 
    for x in $(echo "${runTable}" | awk -F \| '{print $2}' | sort -u)
    do 
    	printf "%15s %25s         %s\n" "JobName" 	"${x}"	
    done    
    printf "%15s %25s         %s\n" "NumJobs" 	"${numJobs}"
    printf "%15s %25s         %s\n" "NumFailedJobs" 	"${failJobs}"
    if [[ ${failJobs} -gt 0 ]]
    then 
        echo "$runTable" | grep -v -e COMPLETED | sed -e "s:^:#:"
    fi
    printf "%15s %25s         %s\n" "CumElapsedTime" 	"${CumElapsedRaw}"         "#sec|min|hours"
    printf "%15s %25s         %s\n" "CumCPUTime"    	"${CumCPUTimeRAW}"         "#sec|min|hours"
    printf "%15s %25s         %s\n" "Partition" 	"${partition}"          "#batch|bigmem|gpu|long|other"
    printf "%15s %25s         %s\n" "ReqMem" 		"${ReqMem}"          	"#Minimum required memory for the job. 'c' at the end - Memory Per CPU, a 'n' - Memory Per Node"
    printf "%15s %25s         %s\n" "NCPUS" 		"${NCPUS}"          	"#Total number of CPUs allocated to the job"
    printf "%15s %25s         %s\n" "AveVMSize" 	"${AveVMSize}"          "#Min|Max|Avg - Average Virtual Memory size of all tasks in job - in Mb"
    printf "%15s %25s         %s\n" "AveRSS" 		"${AveRSS}"          	"#Min|Max|Avg - Average resident set size of all tasks in job - in Mb"
    printf "%15s %25s         %s\n" "NTasks" 		"${NTasks}"          	"#Min|Max|Avg - Total number of tasks in a job or step"
    printf "%15s %25s         %s\n" "AveCPUFreq" 	"${AveCPUFreq}"         "#Min|Max|Avg - Average weighted CPU frequency of all tasks in job - in GHz"
    printf "%15s %25s         %s\n" "ConsumedEnergy" 	"${ConsumedEnergy}"     "#Min|Max|Avg - Total energy consumed by all tasks in job - in joules"
    printf "%15s %25s         %s\n" "AveDiskRead" 	"${AveDiskRead}"        "#Min|Max|Avg - Average number of bytes read by all tasks in job - in Mb"
    printf "%15s %25s         %s\n" "AveDiskWrite" 	"${AveDiskWrite}"       "#Min|Max|Avg - Average number of bytes written by all tasks in job - in Mb"
    printf "%15s %25s         %s\n" "sittingInQueue" 	"$(echo "${Submit} ${End}" | awk '{printf "%d|%1.f|%.1f", $2-$1, ($2-$1)/60, ($2-$1)/3600}')"          "#sec|min|hours  the job is sitting in the queue"    
}

DB=(${RAW_DB} ${FIX_DB} ${COR_DB} ${CONT_DB})
marvelPhases=(fix 10x arrow bionano cont corr filt freebayes hic mask phase purgeHaplotigs scaff10x scrub tour)
marvelJobs=(LAfix arksArks arksLinks arksLongranger arksPrepare arksTigmint arrow bamMerge bamseparate bamtools BNscaffold BNstatistics break10xBreak10x break10xLongrangerAlign break10xPrepare break10xPrepareIntermediate break10xStatistics Catrack createDB createMinimap2RefIndex createSubDir CTanalyze CTstatistics daligner datander DBdust FBbwa FBconsensus FBfastp FBfreebayes FBlongrangerAlign FBmarkDuplicates FBprepareInput FBstatistics forcealign HICbwa HICfilter HiChiglassBwa HiChiglassFilter HiChiglassMatrix HIChiglassPrepare HICmarkduplicates HICmerge HICprepareInput HICsalsaBwa HICsalsaFilter HICsalsaMarkduplicates HICsalsaMerge HICsalsaPrepareInput HICsalsaSalsa HICsalsaStatistics HICstatistics LAcorrect LAfilter LAgap LAmerge LAq LArepeat LAseparate LAstitch LongrangerBcftoolsConsensus LongrangerLongrangerWgs LongrangerPrepareInput LongrangerStatistics marvelStats minimap2 OGbuild OGlayout OGtour paths2rids pbalign prepDB prepInFasta readCovHist repcomp scaff10Xbreak10x scaff10xLongrangerAlign scaff10xprepare scaff10Xprepare scaff10xPrepareIntermediate scaff10xScaff10x scaff10Xscaff10x scaff10xStatistics scaff10Xstatistics statistics TANmask TKcombine TKhomogenize TKmerge tour2fasta)

dbCount=0
while [ "x${DB[dbCount]}" != "x" ]
do
	db=${DB[dbCount]}
	jobCount=0
	while [ "x${marvelJobs[jobCount]}" != "x" ]
    do
		job=${marvelJobs[jobCount]}
		phaseCount=0            	        
        while [ "x${marvelPhases[phaseCount]}" != "x" ]
        do
        	phase=${marvelPhases[phaseCount]}
        	folder=log_${phase}_${job}_${db}
            if [[ -d ${folder} ]] 
            then    
                
                blacklist="" 

                for x in ${folder}/*.out
                do
                    if [[ ! -f $x ]]       ### if folder is empty 
                    then
                        break;
                    fi
                   
                    slurm_id=""
                    if [[ "${x}" =~ ^${folder}/${phase}_([0-9]*)_${db}_([0-9]*).out$ ]]
                    then
                        slurm_id=$(echo ${x%.out} | awk -F _ '{print $NF}')
                    elif [[ "${x}" =~ ^${folder}/${phase}_([0-9]*)_${db}_([0-9]*)_([0-9]*).out$ ]]
                    then
                        slurm_id=$(echo ${x%.out} | awk -F _ '{print $(NF-1)}')
                	elif [[ "${x}" =~ ^${folder}/${phase}_([0-9]*)_${db}_([0-9]*)_([0-9]*)_([0-9]*).out$ ]] # in case: job files were splitted into multiple parts
                    then
                        slurm_id=$(echo ${x%.out} | awk -F _ '{print $(NF-1)}')
                    fi

                    stop=0
                    for y in $blacklist
                    do
                        if [[ "$y" == "$slurm_id" ]]
                        then
                            stop=1
                            break
                        fi
                    done

                    if [[ $stop -eq 1 || -z $slurm_id ]]
                    then
                        continue
                    fi

                    blacklist="$slurm_id ${blacklist}"

                    if [[ ! -s stats/${phase}_${job}_${db}_${slurm_id}.txt ]]
                    then
                        createSlurmStats ${slurm_id} > stats/${phase}_${job}_${db}_${slurm_id}.txt
                        jobname=$(grep -e JobName stats/${phase}_${job}_${db}_${slurm_id}.txt | awk '{print $2}')
                        fext=$(echo ${jobname} | sed -e "s: :_:g")
                		if [[ $(echo ${fext} | wc -m) -gt 0 ]]
                		then
                			ln -s -f -r stats/${phase}_${job}_${db}_${slurm_id}.txt stats/${phase}_${job}_${db}_${slurm_id}_${fext}.lnk	
                		fi
                    fi
                done        
            fi
            phaseCount=$((${phaseCount}+1))
        done
        jobCount=$((${jobCount}+1))
    done
    dbCount=$((${dbCount}+1))
done

### todo
#### 1 update overall stats (if neccessary)
reportVariables=(NumJobs NumFailedJobs CumElapsedTime CumCPUTime sittingInQueue)
# header
header=$(printf "%20s %15s" "marvelPhase" "numSlurmJobs")
reportCount=0   	
while [ "x${reportVariables[reportCount]}" != "x" ]
do
	if [[ ${reportVariables[reportCount]} =~ Num ]]
	then
		header="${header}$(printf " %15s" "${reportVariables[reportCount]}")"
	else
		header="${header}$(printf " %25s" "${reportVariables[reportCount]}[s|m|h]")"
	fi						      			
	reportCount=$(( $reportCount + 1 ))
done        		        

dbCount=0
while [ "x${DB[dbCount]}" != "x" ]
do 
	db=${DB[dbCount]}
    phaseCount=0
    phaseReport=""            	        
    while [ "x${marvelPhases[phaseCount]}" != "x" ]
    do
    	phase=${marvelPhases[phaseCount]}    	    	    	
		jobCount=0		
		while [ "x${marvelJobs[jobCount]}" != "x" ]
        do
        	appliedJobs[${jobCount}]=0
            jobCount=$(( $jobCount + 1 ))
        done		
		jobCount=0
		body=""
		while [ "x${marvelJobs[jobCount]}" != "x" ]
    	do	
    		job=${marvelJobs[jobCount]}
    		found=0
    		# reset reported Vars		 	
	    	reportCount=0   	
	    	while [ "x${reportVariables[reportCount]}" != "x" ]
	        do
	        	appliedReportVariables[${reportCount}]=0
	            reportCount=$(( $reportCount + 1 ))
	        done                
    		    		       	
        	for f in stats/${phase}_${job}_${db}_*.txt
        	do 
        		if [[ -f ${f} ]]
				then
					found=1
            		reportCount=0
            		appliedJobs[${jobCount}]=$((${appliedJobs[jobCount]}+1))        				            		   	
    				while [ "x${reportVariables[reportCount]}" != "x" ]
        			do
        				tmp=$(grep -e ${reportVariables[${reportCount}]} $f | awk '{print $2}' |  awk -F \| '{print $1}')
        				if [[ "x${tmp}" == "x" ]]
        				then
							(>&2 echo "Unable to get variable ${reportVariables[${reportCount}]} from file $f")
            				exit 1
        				fi
						appliedReportVariables[${reportCount}]=$((${appliedReportVariables[reportCount]}+${tmp}))
            			reportCount=$(( $reportCount + 1 ))            			
        			done                
				fi
        	done
   
			# report stats for individual phases 
			if [[ ${found} -eq 1 ]]
			then		   	
		   		body="${body}$(printf "%20s %15s" "${marvelJobs[jobCount]}" "${appliedJobs[jobCount]}")"
				reportCount=0   	
				while [ "x${reportVariables[reportCount]}" != "x" ]
    			do
    				if [[ ${reportVariables[reportCount]} =~ Num ]]
    				then
    					body="${body}$(printf " %15s" "${appliedReportVariables[reportCount]}")"
    				else
    					body="${body}$(printf " %25s" "$(echo ${appliedReportVariables[reportCount]} | awk '{printf "%d|%.1f|%.1f\n", $1, $1/60, $1/3600}')")"
    				fi		        				     				
        			reportCount=$(( $reportCount + 1 ))
    			done   
    			     
    			body="${body}\n"        							    	
			fi
       		jobCount=$(( $jobCount + 1 ))
    	done 
		if [[ -n ${body} ]]
		then
			printf "%b\n" "${header}" > stats/all_${phase}_${db}.txt
			printf "%b" "${body}" >> stats/all_${phase}_${db}.txt
			### sum up all steps of current phase 
			tmpPhase="$(printf "%20s %15s" "ALL_${phase}_PHASES" "$(printf "%b" "${body}" | awk '{s+=$2} END {print s}')")"
			reportCount=0   	
			while [ "x${reportVariables[reportCount]}" != "x" ]
			do
				if [[ ${reportVariables[reportCount]} =~ Num ]]
				then
					tmpPhase="${tmpPhase}$(printf " %15s" "$(printf "%b" "${body}" | awk -v off=${reportCount} '{s+=$(3+off)} END {print s}')")"					
				else
					tmpPhase="${tmpPhase}$(printf " %25s" "$(printf "%b" "${body}" | awk -v off=${reportCount} '{print $(3+off)}' | awk -F \| '{s+=$1; m+=$2; h+=$3} END {printf "%d|%.1f|%.1f\n", s, m, h}')")"					
				fi		        				     				
				reportCount=$(( $reportCount + 1 ))
			done   
			tmpPhase="${tmpPhase}\n"
			phaseReport="${phaseReport}${tmpPhase}"
			printf "%b" "${tmpPhase}" >> stats/all_${phase}_${db}.txt		
		fi 
	    phaseCount=$(( $phaseCount + 1 ))   	                    
    done
	printf "%b\n" "${header}" > stats/all_phases_${db}.txt
	printf "%b" "${phaseReport}" >> stats/all_phases_${db}.txt	
	### sum up all phases
	tmpPhase=$(printf "%20s %15s" "ALL_PHASES" "$(printf "%b" "${phaseReport}" | awk '{s+=$2} END {print s}')")
	reportCount=0   	
	while [ "x${reportVariables[reportCount]}" != "x" ]
	do
		if [[ ${reportVariables[reportCount]} =~ Num ]]
		then
			tmpPhase="${tmpPhase}$(printf " %15s" "$(printf "%b" "${phaseReport}" | awk -v off=${reportCount} '{s+=$(3+off)} END {print s}')")"					
		else
			tmpPhase="${tmpPhase}$(printf " %25s" "$(printf "%b" "${phaseReport}" | awk -v off=${reportCount} '{print $(3+off)}' | awk -F \| '{s+=$1; m+=$2; h+=$3} END {printf "%d|%.1f|%.1f\n", s, m, h}')")"					
		fi		        				     				
		reportCount=$(( $reportCount + 1 ))
	done   
	tmpPhase="${tmpPhase}\n"
	printf "%b" "${tmpPhase}" >> stats/all_phases_${db}.txt     
    dbCount=$(( $dbCount + 1 ))
done