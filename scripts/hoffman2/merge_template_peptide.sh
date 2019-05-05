#!/bin/bash
#$ -cwd
#$ -o ../output/merge_temp_pep.txt
#$ -j y
#$ -l h_data=10000M,h_rt=3:00:00
#$ -pe shared 8
#$ -t 2-5:1

. /u/local/Modules/default/init/modules.sh
module load python/anaconda3

source activate rdkit

python ./merge_template_peptide.py -pin length${SGE_TASK_ID}_all.txt