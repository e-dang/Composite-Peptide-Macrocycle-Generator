from argparse import ArgumentParser

import macrocycles.ranges as ranges
import macrocycles.project_io as project_io
import macrocycles.runners as runners

parser = ArgumentParser(description='Parses arguments to determine what chunk of template_peptides to work on.')
parser.add_argument('--peptide_len', type=int, required=True, help='The length of peptide to work on.')
parser.add_argument('--num_jobs', type=int, default=1, help='The total number of jobs requested for this part. Can '
                    'be thought of as how many chunks the whole set of macrocycles should be split up into.')
parser.add_argument('--num', type=int, default=1, help='The job number assigned to this process, i.e. what chunk '
                    'should this process work on.')
parser.add_argument('--num_conformers', type=int, required=True, help='The number of conformers to generate plan for.')
parser.add_argument('--ebejer', action='store_true', help='Determines whether to use Ebejer method or not.')

args = parser.parse_args()

data_chunk = ranges.ContinuousDataChunk(args.num_conformers, args.num_jobs, args.num)

if args.ebejer:
    runners.run_ebejer(peptide_length=args.peptide_len, data_chunk=data_chunk, job_num=args.num)
else:
    runners.run_conformers(peptide_length=args.peptide_len, data_chunk=data_chunk, job_num=args.num)
