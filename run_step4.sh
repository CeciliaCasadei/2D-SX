#!/bin/bash


#SBATCH --cpus 10
#SBATCH --job-name=step4


module load psi-python27

python runSetup_4_mergeVsModel.py


