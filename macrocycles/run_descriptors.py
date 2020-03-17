import math
import os
from argparse import ArgumentParser

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


parser = ArgumentParser(description='Parses arguments to determine which descriptor to calculate')
parser.add_argument('--mw', action='store_true', help='Calculate molecular weight descriptor')
parser.add_argument('--rb', action='store_true', help='Calculate rotatable bonds descriptor')
parser.add_argument('--tpsa', action='store_true', help='Calculate TPSA descriptor')
parser.add_argument('--peptide_len', type=int, required=True, help='The length of peptide in the macrocycles.')
parser.add_argument('--num_jobs', type=int, default=1, help='The total number of jobs requested for this part. Can '
                    'be thought of as how many chunks the whole set of macrocycles should be split up into.')
parser.add_argument('--num', type=int, default=1, help='The job number assigned to this process, i.e. what chunk '
                    'should this process work on.')
parser.add_argument('--macrocycle_output', type=str, help='The path to the file containing the output of '
                    'macrocycle generation step.')
parser.add_argument('--num_macrocycles', type=int, help='The number of macrocycles to calculate descriptors for.')

args = parser.parse_args()

if sum([args.mw, args.rb, args.tpsa]) != 1:
    parser.error('Must specify ONE descriptor!')

if args.num_macrocycles is None:
    num_macrocycles = get_num_macrocycles(args.macrocycle_output, args.peptide_len)
else:
    num_macrocycles = args.num_macrocycles

# calculate the chunk to work on
chunk_size = math.ceil(num_macrocycles / args.num_jobs)
start = chunk_size * (args.num - 1)
end = chunk_size * args.num
if end > num_macrocycles:
    end = num_macrocycles


elif args.mw:
    runners.run_mw_descriptor(peptide_length=args.peptide_len, start=start, end=end, job_num=args.num)
elif args.rb:
    runners.run_rb_descriptor(peptide_length=args.peptide_len, start=start, end=end, job_num=args.num)
else:
    runners.run_tpsa_descriptor(peptide_length=args.peptide_len, start=start, end=end, job_num=args.num)
