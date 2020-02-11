from argparse import ArgumentParser

import macrocycles.runners as runners
import os
import math


def get_num_macrocycles(filepath):
    count = 0
    for file_name in os.listdir(filepath):
        if 'part3' in file_name:
            with open(os.path.join(filepath, file_name), encoding='windows-1252') as file:
                for line in file.readlines():
                    line = line.split(' ')
                    count += int(line[-1].strip('\n'))

    return count


parser = ArgumentParser(description='Parses arguments to determine what chunk of template_peptides to work on.')
parser.add_argument('--peptide_len', type=int, required=True, help='The length of peptide to work on.')
parser.add_argument('--num_jobs', type=int, default=1, help='The total number of jobs requested for this part. Can '
                    'be thought of as how many chunks the whole set of macrocycles should be split up into.')
parser.add_argument('--num', type=int, default=1, help='The job number assigned to this process, i.e. what chunk '
                    'should this process work on.')
parser.add_argument('--macrocycle_output', type=int, help='The path to the stdout file during '
                    'macrocycle generation.')
parser.add_argument('--num_macrocycles', type=int, help='The number of macrocycles generated in last '
                    'round.')

args = parser.parse_args()

# get total number of template_peptides with specified peptide length
if args.num_macrocycles:
    num_macrocycles = args.num_macrocycles
elif args.macrocycle_output:
    num_macrocycles = get_num_macrocycles(args.macrocycle_output)
else:
    parser.error('Must specify either option --macrocycle_output or --num_macrocycles.')

# calculate the chunk to work on
chunk_size = math.ceil(num_macrocycles / args.num_jobs)
start = chunk_size * (args.num - 1)
end = chunk_size * args.num
if end > num_macrocycles:
    end = num_macrocycles

runners.run_conformers(peptide_length=args.peptide_len, start=start, end=end, job_num=args.num)
