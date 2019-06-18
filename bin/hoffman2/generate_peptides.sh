#!/bin/bash
#$ -cwd
#$ -o ../output/gen_peptides.txt
#$ -j y
#$ -l h_data=12000M,h_rt=10:00:00
#$ -pe shared 4
#$ -t 1-1000:1

. /u/local/Modules/default/init/modules.sh
module load python/anaconda3

source activate rdkit

python ./generate_peptides.py 2 --num_jobs 1000 --job_num $SGE_TASK_ID --fp_out /u/scratch/e/ericdang/peptides/