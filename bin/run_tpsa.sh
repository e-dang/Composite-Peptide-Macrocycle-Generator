#!/bin/bash
#$ -cwd
#$ -o ../output/tpsa.txt
#$ -j y
#$ -l h_data=1000M,h_rt=10:00:00,h_vmem=8000M,highp
#$ -pe shared 8
#$ -t 3-4:1

. /u/local/Modules/default/init/modules.sh
module load python/anaconda3

source activate rdkit

python ../macrocycles/run_descriptors.py --tpsa --peptide_len ${SGE_TASK_ID}