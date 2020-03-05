from argparse import ArgumentParser
from time import time
import os

import macrocycles.runners as runners


def get_num_macrocycles(filepath, peptide_len):
    count = 0
    for file_name in os.listdir(filepath):
        if 'part3_' + str(peptide_len) in file_name:
            with open(os.path.join(filepath, file_name), encoding='windows-1252') as file:
                for line in file.readlines():
                    line = line.split(' ')
                    count += int(line[-1].strip('\n'))
    return count


parser = ArgumentParser(description='Parses arguments to determine what chunk of template_peptides to work on.')
parser.add_argument('--peptide_len', type=int, required=True, help='The length of peptide to work on.')
parser.add_argument('--num_conformers', type=int, required=True, help='The number of conformers to generate plan for.')
parser.add_argument('--macrocycle_output', type=str, help='The path to the file containing the output of '
                    'macrocycle generation step.')

args = parser.parse_args()

if args.macrocycle_output:
    num_macrocycles = get_num_macrocycles(args.macrocycle_output, args.peptide_len)
else:
    num_macrocycles = None

start = time()
runners.generate_conformer_plan(args.peptide_len, args.num_conformers, num_macrocycles)
print(
    f'Completed conformer plan with peptide length {args.peptide_len} for {args.num_conformers} conformers in {time() - start}')
