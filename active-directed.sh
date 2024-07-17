#!/bin/bash
source configs/run.env
shopt -s globstar
lockdir=/tmp/mcsp_ll.lock

write_output () {
  if mkdir -- "$lockdir" 2> /dev/null
  then
    echo ${@:1:$#-1} >> ${@: -1}
    rm -rf $lockdir
    return 0
  else
    return 1
  fi
}

task () {
  out1=$(./bin/run.o min_max ./${1} ./${2} -l -q -i -t ${TIMEOUT})
  out1="${1}, ${2}, ${out1}"

  available=1
  while [ "$available" -eq 1 ]
  do
    write_output ${out1} ${3}
    available=$?
  done
  return
}

# Files to process
files=("biochemicalReactions_diff")

# Maximum number of parallel tasks

# Function to wait for available slot
wait_for_slot() {
  while [ $(jobs -r | wc -l) -ge $N_PARALLEL ]; do
    sleep 1
  done
}

# Iterate over the files
for f in ${files[@]}; do
  file_name=indicators/${f}_data.txt
  results_file_name=results/${f}_${MCSP_HEURISTIC}_${TIMEOUT}.txt

  let j=0
  if [ -e $results_file_name ]; then
    while IFS= read -r line; do
      j=$((j + 1))
    done <"$results_file_name"
  fi

  while IFS=',' read -r file1 file2; do
    if [[ $j -ne 0 ]]; then
      j=$((j - 1))
      continue
    fi

    wait_for_slot
    task $file1 $file2 $results_file_name &

  done <"$file_name"
done

# Wait for all tasks to complete
wait
