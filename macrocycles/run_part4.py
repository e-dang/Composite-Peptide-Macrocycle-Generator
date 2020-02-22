from argparse import ArgumentParser
from time import time

import macrocycles.runners as runners

parser = ArgumentParser(description='Parses arguments to determine what chunk of template_peptides to work on.')
parser.add_argument('--peptide_len', type=int, required=True, help='The length of peptide to work on.')
parser.add_argument('--num_macrocycles', type=int, help='The number of macrocycles generated in last '
                    'round.')

args = parser.parse_args()

start = time()
runners.generate_conformer_plan(args.peptide_len, args.num_macrocycles)
print(
    f'Completed conformer plan with peptide length {args.peptide_len} for {args.num_macrocycles} macrocycles in {time() - start}')
