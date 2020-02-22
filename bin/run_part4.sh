#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_data=20000M,h_rt=05:00:00
#$ -t 3-5:1
#$ -o ../output/part4_${SGE_TASK_ID}.txt


. /u/local/Modules/default/init/modules.sh
module load python/anaconda3

source activate rdkit

python ../macrocycles/run_part4.py --peptide_len ${SGE_TASK_ID} --num_macrocycles 333333