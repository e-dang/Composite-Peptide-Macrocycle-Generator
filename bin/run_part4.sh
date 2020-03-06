#!/bin/bash
#$ -cwd
#$ -o ../output/part4.txt
#$ -j y
#$ -l h_data=1000M,h_rt=00:10:00,h_vmem=1000M
#$ -t 3-5:1


. /u/local/Modules/default/init/modules.sh
module load python/anaconda3

source activate rdkit

if [ $SGE_TASK_ID -eq 3 ]
then
    NUM_CONFS=333334
elif [ $SGE_TASK_ID -eq 4 ]
then
    NUM_CONFS=333333
else
    NUM_CONFS=333333
fi

python ../macrocycles/run_part4.py --peptide_len ${SGE_TASK_ID} --num_conformers $NUM_CONFS --macrocycle_output ../output/ > ../output/part4_${SGE_TASK_ID}.txt