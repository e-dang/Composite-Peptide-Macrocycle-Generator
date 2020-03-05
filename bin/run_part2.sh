#!/bin/bash
#$ -cwd
#$ -o ../output/part2.txt
#$ -j y
#$ -l h_data=2000M,h_rt=18:00:00,h_vmem=16000M,highp
#$ -pe shared 8
#$ -t 3-5:1

. /u/local/Modules/default/init/modules.sh
module load python/anaconda3

source activate rdkit

if [ $SGE_TASK_ID -eq 3 ]
then
    NUM_PEPTIDES=1500000
elif [ $SGE_TASK_ID -eq 4 ]
then
    NUM_PEPTIDES=3500000
else
    NUM_PEPTIDES=1500000
fi

python ../macrocycles/run_part2.py --peptide_len $SGE_TASK_ID --num_peptides $NUM_PEPTIDES > ../output/part2_${SGE_TASK_ID}.txt