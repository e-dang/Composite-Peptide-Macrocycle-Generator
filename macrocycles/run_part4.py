from argparse import ArgumentParser

import macrocycles.runners as runners

parser = ArgumentParser(description='Parses arguments to determine what chunk of template_peptides to work on.')
parser.add_argument('--peptide_len', type=int, required=True, help='The length of peptide to work on.')

args = parser.parse_args()

runners.run_conformers(peptide_length=args.peptide_len)
