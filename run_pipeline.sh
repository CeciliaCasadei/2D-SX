#!/bin/bash


jid1=$(sbatch -p long --job-name=step1 --parsable run_step1_3.sh runSetup_1.py)

# multiple jobs can depend on a single job
jid2=$(sbatch -p long --job-name=step2 --parsable --dependency=afterok:$jid1 run_step2.sh)

jid3=$(sbatch -p long --job-name=step3 --parsable --dependency=afterok:$jid2 run_step1_3.sh runSetup_3_modelling.py)

jid4=$(sbatch -p long  --job-name=step4 --parsable --dependency=afterok:$jid3 run_step4.sh)
