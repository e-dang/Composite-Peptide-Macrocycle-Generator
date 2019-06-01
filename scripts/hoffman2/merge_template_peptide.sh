#!/bin/bash
#$ -cwd
#$ -o ../output/merge_temp_pep.txt
#$ -j y
#$ -l h_data=12000M,h_rt=5:00:00
#$ -t 1-1000:1

. /u/local/Modules/default/init/modules.sh
module load python/anaconda3

source activate rdkit

python ./merge_template_peptide.py --fp_pep /u/scratch/e/ericdang/peptides/ --f_pep length2_all_${SGE_TASK_ID}*.json --fp_out /u/scratch/e/ericdang/template_peptides/