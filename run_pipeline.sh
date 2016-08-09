#!/bin/bash


jid1=$(sbatch --job-name=step1 --parsable run_step1_3.sh runSetup_1.py)

# multiple jobs can depend on a single job
#jid2=$(sbatch --job-name=step2 --parsable --dependency=afterok:$jid1 run_step2.sh)

#jid3=$(sbatch --job-name=step3 --parsable --dependency=afterok:$jid2 run_step1_3.sh runSetup_3_modelling.py)
