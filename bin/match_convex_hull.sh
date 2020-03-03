#!/bin/bash
#$ -cwd
#$ -o ../output/convex_hull.txt
#$ -j y
#$ -l h_data=1000M,h_rt=01:00:00,highp

. /u/local/Modules/default/init/modules.sh
module load python/anaconda3

source activate rdkit

python ../macrocycles/match_convex_hull.py