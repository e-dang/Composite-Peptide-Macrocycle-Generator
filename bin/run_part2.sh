#!/bin/bash
#$ -cwd
#$ -o ../output/part2.txt
#$ -j y
#$ -l h_data=1000M,h_rt=05:00:00
#$ -pe shared 8
#$ -t 3-5:1

. /u/local/Modules/default/init/modules.sh
module load python/anaconda3

source activate rdkit

python ../macrocycles/run_part2.py --peptide_len $SGE_TASK_ID > ../output/part2_${SGE_TASK_ID}.txt