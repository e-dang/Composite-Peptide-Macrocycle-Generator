#!/bin/bash
#$ -cwd
#$ -o ../output/part5_5.txt
#$ -j y
#$ -l h_data=3000M,h_rt=05:00:00
#$ -pe shared 8
#$ -t 1-1500:1

. /u/local/Modules/default/init/modules.sh
module load python/anaconda3
module load openbabel/2.4.1

source activate rdkit

python ../macrocycles/run_part5.py --peptide_len 5 --num_jobs 1500 --num $SGE_TASK_ID > ../output/part5_5_${SGE_TASK_ID}.txt