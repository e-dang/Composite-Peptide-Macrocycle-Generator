#!/bin/bash

qsub -N part1 run_part1.sh
qsub -N part2 -hold_jid part1 run_part2.sh
qsub -N part3_3 -hold_jid part2 run_part3_3.sh
qsub -N part3_4 -hold_jid part2 run_part3_4.sh
qsub -N part3_5 -hold_jid part2 run_part3_5.sh