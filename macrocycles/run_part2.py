from argparse import ArgumentParser

import macrocycles.runners as runners

NUM_PEPTIDES = 3000

parser = ArgumentParser(description='Parse arguments to determine the length of peptides to generate.')

parser.add_argument('--peptide_len', type=int, required=True, help='The length of the peptides to generate.')

args = parser.parse_args()

# create peptide plan
runners.generate_peptide_plan(args.peptide_len, NUM_PEPTIDES)

# create peptides
runners.run_peptides(peptide_length=args.peptide_len)

# create template_peptides
runners.run_template_peptides(peptide_length=args.peptide_len)
