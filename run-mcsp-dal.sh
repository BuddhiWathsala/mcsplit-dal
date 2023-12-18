#!/bin/bash
source configs/run.env
shopt -s globstar
lockdir=/tmp/mcsp_ll.lock
write_output () {
  if mkdir -- "$lockdir" 2> /dev/null
  then
    # echo ${@: -1}
    echo ${@:1:$#-1} >> ${@: -1}
    rm -rf $lockdir
    return 0
  else
    return 1
  fi
}

task () {
  # echo ./bin/run.o ./${1} ./${2} -l -t ${TIMEOUT} -i ${INIT_COUNT}
  # out1=$(/Users/buddhi/Documents/Work/PhD/Implementation/mcsplit-plus/bin/run.o ./${1} ./${2} -l -q -t ${TIMEOUT} -i ${INIT_COUNT} -m ${MAX_INNER_THRESHOLD})
  out1=$(./bin/run.o min_max ./${1} ./${2}  -c -l -q -t ${TIMEOUT})
  # echo ${MAX_INNER_THRESHOLD}
  out1="${1}, ${2}, ${out1}"
#   out1=${out1//[,]/' '} 
  
  available=1
  while [ "$available" -eq 1 ]
  do
    write_output ${out1} ${3}
    available=$?
  done
  return
}

# files=("largerLV" "si" "scalefree" "scalefree" "meshes-CVIU11" "LV" "images-CVIU11" "images-PR15")
#files=("completedLV")
files=("LV" "largerLV" "meshes-CVIU11"  "images-CVIU11" "si" "scalefree" "phase")
for f in ${files[@]}; do
    file_name=indicators/${f}_data.txt
    results_file_name=results/${f}_${MCSP_HEURISTIC}_${TIMEOUT}.txt
    # rm -f $results_file_name
    let i=0
    let j=0
    if [ -e $results_file_name ]
    then
      while IFS= read -r line
      do 
          j=$((j + 1))
      done <"$results_file_name"
    fi
    
    while IFS=',' read -r file1 file2
    do
        if [[ $j -ne 0 ]]
        then
          j=$((j - 1))
          continue
        fi
        # echo $file1 $file2 $results_file_name
        task $file1 $file2 $results_file_name &
        if [[ ($i -ne 0) && ($i%$N_PARALLEL -eq 0) ]]
        then
            wait
        fi
        i=$((i + 1))
    done <"$file_name"
done



