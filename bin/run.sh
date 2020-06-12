#!/bin/bash

ALLOC_HOURS=1
ALLOC_MINUTES=0
ALLOC_SECONDS=0
ALLOC_TIME=$((${ALLOC_HOURS} * 60 * 60 + ${ALLOC_MINUTES} * 60 + ${ALLOC_SECONDS}))
SECONDS=0


echo "Initializing..." &&
./cpmg.sh initialize &&
echo "Inactivating extra sidechains..." &&
./cpmg.sh deactivate sidechain -f /Users/ericdang/Documents/Research/Macrocycles/data/chemdraw/conservative_sidechains.sdf -i
echo "Generating sidechains..." &&
./cpmg.sh generate -p multi -t $((${ALLOC_TIME} - ${SECONDS})) sidechain &&
echo "Generating monomers..." &&
./cpmg.sh generate -p multi -t $((${ALLOC_TIME} - ${SECONDS})) monomer &&
echo "Generating inter reactions..." &&
./cpmg.sh generate -t $((${ALLOC_TIME} - ${SECONDS})) inter_reaction &&
echo "Generating intra reactions..." &&
./cpmg.sh generate -t $((${ALLOC_TIME} - ${SECONDS})) intra_reaction &&
echo "Generating peptide plan for peptides of length 3..." &&
./cpmg.sh generate -t $((${ALLOC_TIME} - ${SECONDS})) peptide_plan -l 3 -n 100 &&
echo "Generating peptide plan for peptides of length 4..." &&
./cpmg.sh generate -t $((${ALLOC_TIME} - ${SECONDS})) peptide_plan -l 4 -n 100 &&
echo "Generating peptide plan for peptides of length 5..." &&
./cpmg.sh generate -t $((${ALLOC_TIME} - ${SECONDS})) peptide_plan -l 5 -n 100 &&
echo "Generating peptides of length 4..." &&
./cpmg.sh generate -p multi -t $((${ALLOC_TIME} - ${SECONDS})) peptide -l 3 &&
echo "Generating peptides of length 4..." &&
./cpmg.sh generate -p multi -t $((${ALLOC_TIME} - ${SECONDS})) peptide -l 4 &&
echo "Generating peptides of length 5..." &&
./cpmg.sh generate -p multi -t $((${ALLOC_TIME} - ${SECONDS})) peptide -l 5 &&
echo "Generating template peptides of length 3..." &&
./cpmg.sh generate -p multi -t $((${ALLOC_TIME} - ${SECONDS})) template_peptide -l 3 &&
echo "Generating template peptides of length 4..." &&
./cpmg.sh generate -p multi -t $((${ALLOC_TIME} - ${SECONDS})) template_peptide -l 4 &&
echo "Generating template peptides of length 5..." &&
./cpmg.sh generate -p multi -t $((${ALLOC_TIME} - ${SECONDS})) template_peptide -l 5 &&
echo "Generating macrocycles of length 3..." &&
./cpmg.sh generate -p multi -t $((${ALLOC_TIME} - ${SECONDS})) macrocycle -l 3
echo "Generating macrocycles of length 4..." &&
./cpmg.sh generate -p multi -t $((${ALLOC_TIME} - ${SECONDS})) macrocycle -l 4
echo "Generating macrocycles of length 5..." &&
./cpmg.sh generate -p multi -t $((${ALLOC_TIME} - ${SECONDS})) macrocycle -l 5
echo "Generating conformers..." &&
./cpmg.sh generate -p multi -t $((${ALLOC_TIME} - ${SECONDS})) conformer -l 3
echo "Generating conformers..." &&
./cpmg.sh generate -p multi -t $((${ALLOC_TIME} - ${SECONDS})) conformer -l 4
echo "Generating conformers..." &&
./cpmg.sh generate -p multi -t $((${ALLOC_TIME} - ${SECONDS})) conformer -l 5