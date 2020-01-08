#!/bin/bash
#$ -cwd
#$ -o ../output/part3_5.txt
#$ -j y
#$ -l h_data=750M,h_rt=02:00:00
#$ -pe shared 8
#$ -t 1-500:1

. /u/local/Modules/default/init/modules.sh
module load python/anaconda3

source activate rdkit

python ../macrocycles/run_part3.py --peptide_len 5 --num_jobs 500 --num $SGE_TASK_ID --tp_out ../output/part2_5.txt > ../output/part3_5_${SGE_TASK_ID}.txt