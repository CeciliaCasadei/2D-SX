#!/bin/bash

#SBATCH --cpus 10
#SBATCH --job-name=step2

cmd=echo

module load psi-python27

$cmd python runSetup_2_merging.py




