#!/bin/bash
JOB1="gen_pep"

qsub -N $JOB1 generate_peptides.sh

qsub -hold_jid $JOB1 merge_template_peptide.sh