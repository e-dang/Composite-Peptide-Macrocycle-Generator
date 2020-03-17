#!/bin/bash
#$ -cwd
#$ -N p5_4
#$ -o ../output/part5_4.txt
#$ -j y
#$ -l h_data=1200M,h_rt=05:00:00,h_vmem=9600M
#$ -pe shared 8
#$ -t 1-1500:1

. /u/local/Modules/default/init/modules.sh
module load python/anaconda3
module load openbabel/2.4.1

source activate rdkit

python ../macrocycles/run_part5.py --peptide_len 4 --num_jobs 1500 --num $SGE_TASK_ID > ../output/part5_4_${SGE_TASK_ID}.txt