#!/bin/bash
#$ -cwd
#$ -o ../output/part1.txt
#$ -j y
#$ -l h_data=1000M,h_rt=00:10:00

. /u/local/Modules/default/init/modules.sh
module load python/anaconda3

source activate rdkit

python ../macrocycles/run_part1.py