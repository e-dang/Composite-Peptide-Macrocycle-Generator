#!/bin/bash
#$ -cwd
#$ -o ../output/part2.txt
#$ -j y
#$ -l h_data=5000M,h_rt=05:00:00

. /u/local/Modules/default/init/modules.sh
module load python/anaconda3

source activate rdkit

python ../macrocycles/run_part2.py