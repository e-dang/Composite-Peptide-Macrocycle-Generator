#!/bin/bash
#$ -cwd
#$ -o ../output/rb.txt
#$ -j y
#$ -l h_data=1000M,h_rt=10:00:00,h_vmem=8000M
#$ -pe shared* 8

. /u/local/Modules/default/init/modules.sh
module load python/anaconda3

source activate rdkit

python ../macrocycles/run_descriptors.py --rb