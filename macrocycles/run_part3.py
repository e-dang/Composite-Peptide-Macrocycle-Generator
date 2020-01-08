import math
from argparse import ArgumentParser

import runners
from utils import suppress_stdout_stderr

parser = ArgumentParser(description='Parses arguments to determine what chunk of template_peptides to work on.')
parser.add_argument('--peptide_len', type=int, required=True, help='The length of peptide to work on.')
parser.add_argument('--num_jobs', type=int, required=True, help='The total number of jobs requested for this part. Can '
                    'be thought of as how many chunks the whole set of template_peptides should be split up into.')
parser.add_argument('--num', type=int, required=True, help='The job number assigned to this process, i.e. what chunk '
                    'should this process work on.')
parser.add_argument('--tp_output', type=str, required=True, help='The path to the file containing the total number of '
                    'template_peptides with a given peptide length.')

args = parser.parse_args()

# get total number of template_peptides with specified peptide length
with open(args.tp_output, 'r') as file:
    for line in file.readlines():
        if 'Template Peptides' in line:
            split_line = line.split(' ')
            if int(split_line[2].strip(':')) == args.peptide_len:
                num_template_peptides = int(split_line[-1].strip('\n'))
                break

# calculate the chunk to work on
chunk_size = math.ceil(num_template_peptides / args.num_jobs)
start = chunk_size * (args.num - 1)
end = chunk_size * args.num
if end > num_template_peptides:
    end = num_template_peptides

with suppress_stdout_stderr():
    print_str = runners.run_macrocycles(peptide_length=args.peptide_len, start=start, end=end)

print(*print_str)
