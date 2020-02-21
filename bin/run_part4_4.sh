#!/bin/bash
#$ -cwd
#$ -o ../output/part4_3.txt
#$ -j y
#$ -l h_data=1000M,h_rt=05:00:00
#$ -pe shared 8
#$ -t 1-1500:1

. /u/local/Modules/default/init/modules.sh
module load python/anaconda3

source activate rdkit

python ../macrocycles/run_part4.py --peptide_len 4 --num_jobs 1500 --num $SGE_TASK_ID --num_macrocycles 333333 > ../output/part4_3_${SGE_TASK_ID}.txt