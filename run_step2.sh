#!/bin/bash

#SBATCH --cpus 10
#SBATCH --job-name=step2



module load psi-python27

python runSetup_2_merging.py




