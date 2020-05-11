from argparse import ArgumentParser
import cpmg.generators as g
import cpmg.data_handlers as h
import cpmg.parallelizers as p


def execute(operation_string, parallelizer_string, *args):
    generator = g.create_generator_from_string(operation_string, *args)
    data_handler = h.create_handler_from_string(operation_string)
    parallelizer = p.create_parallelizer_from_string(parallelizer_string)
    return parallelizer.execute(generator, data_handler)


def execute_general(args):
    return execute(args.operation, args.parallelism)


def execute_peptide_plan(args):
    return execute(args.operation, args.parallelism, args.length, args.num)


def execute_with_peptide_length(args):
    return execute(args.operation, args.parallelism, args.length)


class GenerateArgParser:
    def __init__(self):
        parser = ArgumentParser()
        subparsers = parser.add_subparsers(dest='operation')
        parser.add_argument('-p', '--parallelism', choices=p.get_all_parallelizer_strings(), nargs='?',
                            const=p.SingleProcess.STRING, default=p.SingleProcess.STRING,
                            help='Selects which level of parallelism to execute the molecule generation with')

        parent_parser = ArgumentParser(add_help=False)
        parent_parser.add_argument('-l', '--length', '--peptide_length', type=int, choices=[3, 4, 5],
                                         required=True,
                                         help='The number of monomers to assemble into a peptide.')

        sidechain_parser = subparsers.add_parser('sidechain')
        sidechain_parser.set_defaults(func=execute_general)

        monomer_parser = subparsers.add_parser('monomer')
        monomer_parser.set_defaults(func=execute_general)

        inter_reaction_parser = subparsers.add_parser('inter_reaction')
        inter_reaction_parser.set_defaults(func=execute_general)

        intra_reaction_parser = subparsers.add_parser('intra_reaction')
        intra_reaction_parser.set_defaults(func=execute_general)

        peptide_plan_parser = subparsers.add_parser('peptide_plan', parents=[parent_parser])
        peptide_plan_parser.add_argument('-n', '--num', '--num_peptides', type=int, required=True,
                                         help='The number of peptides to generate.')
        peptide_plan_parser.set_defaults(func=execute_peptide_plan)

        peptide_parser = subparsers.add_parser('peptide', parents=[parent_parser])
        peptide_parser.set_defaults(func=execute_peptide_plan)

        args = parser.parse_args()
        self.return_val = args.func(args)


if __name__ == "__main__":
    generate_parser = GenerateArgParser()
