from argparse import ArgumentParser

import macrocycles.runners as runners
from macrocycles.utils import suppress_stdout_stderr
import macrocycles.ranges as ranges


def get_num_tp_hybrids(filepath):
    with open(filepath, 'r') as file:
        for line in file.readlines():
            if 'Template Peptides' in line:
                split_line = line.split(' ')
                if int(split_line[2].strip(':')) == args.peptide_len:
                    return int(split_line[-1].strip('\n'))


parser = ArgumentParser(description='Parses arguments to determine what chunk of template_peptides to work on.')
parser.add_argument('--peptide_len', type=int, required=True, help='The length of peptide to work on.')
parser.add_argument('--num_jobs', type=int, default=1, help='The total number of jobs requested for this part. Can '
                    'be thought of as how many chunks the whole set of template_peptides should be split up into.')
parser.add_argument('--num', type=int, default=1, help='The job number assigned to this process, i.e. what chunk '
                    'should this process work on.')
parser.add_argument('--tp_output', type=str, help='The path to the file containing the total number of '
                    'template_peptides with a given peptide length.')
parser.add_argument('--num_tp_hybrids', type=int, help='The number of generated template peptides with a given peptide '
                    'length')

args = parser.parse_args()

# get total number of template_peptides with specified peptide length
if args.num_tp_hybrids:
    num_template_peptides = args.num_tp_hybrids
elif args.tp_output:
    num_template_peptides = get_num_tp_hybrids(args.tp_output)
else:
    parser.error('Must specify either option --tp_output or --num_tp_hybrids.')

data_chunk = ranges.ContinuousDataChunk(num_template_peptides, args.num_jobs, args.num)

with suppress_stdout_stderr():
    print_str = runners.run_macrocycles(peptide_length=args.peptide_len, data_chunk=data_chunk, job_num=args.num)

print(*print_str)
