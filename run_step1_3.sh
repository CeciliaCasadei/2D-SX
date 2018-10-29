#!/bin/bash

#SBATCH --cpus 7
#SBATCH --job-name=step1   # ?????


if [ $# != 1 ]; then
    echo USAGE: $0 script_name
    exit
fi

script=$1

module load psi-python27

if [ ! -d log ]; then
    mkdir log
fi

for run in `seq 195 201`; do
    python $script $run > log/${SLURM_JOB_NAME}_${run}.log &
    echo $run
    # sleep 60 &
done

wait