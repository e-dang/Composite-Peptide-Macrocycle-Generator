from argparse import ArgumentParser

import macrocycles.runners as runners

parser = ArgumentParser(description='Parses arguments to determine which descriptor to calculate')
parser.add_argument('--mw', action='store_true', help='Calculate molecular weight descriptor')
parser.add_argument('--rb', action='store_true', help='Calculate rotatable bonds descriptor')
parser.add_argument('--tpsa', action='store_true', help='Calculate TPSA descriptor')

args = parser.parse_args()

if sum([args.mw, args.rb, args.tpsa]) != 1:
    parser.error('Must specify ONE descriptor!')
elif args.mw:
    runners.run_mw_descriptor(peptide_length=None, job_num=None)
elif args.rb:
    runners.run_rb_descriptor(peptide_length=None, job_num=None)
else:
    runners.run_tpsa_descriptor(peptide_length=None, job_num=None)
