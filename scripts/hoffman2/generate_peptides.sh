#!/bin/bash
#$ -cwd
#$ -o ../output/gen_peptides.txt
#$ -j y
#$ -l h_data=12000M,h_rt=5:00:00
#$ -pe shared 8
#$ -t 2-5:1

. /u/local/Modules/default/init/modules.sh
module load python/anaconda3

source activate rdkit

python ./generate_peptides.py $SGE_TASK_ID