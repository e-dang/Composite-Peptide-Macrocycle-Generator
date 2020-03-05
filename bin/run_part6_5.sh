#!/bin/bash
#$ -cwd
#$ -N p6_5
#$ -o ../output/part6_5.txt
#$ -j y
#$ -l h_data=1200M,h_rt=05:00:00,h_vmem=9600M
#$ -pe shared 8
#$ -t 1-1500:1

. /u/local/Modules/default/init/modules.sh
module load python/anaconda3

source activate rdkit

python ../macrocycles/run_part5.py --peptide_len 5 --num_jobs 1500 --num $SGE_TASK_ID --ebejer > ../output/part6_5_${SGE_TASK_ID}.txt