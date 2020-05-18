from argparse import ArgumentParser
import cpmg.parallelism as p
import cpmg.orchestrator as o
import cpmg.generators as g


def execute(command_line_args):
    p.Parallelism.set_level(command_line_args.parallelism)
    params = o.ExecutionParameters(vars(command_line_args))
    orchestrator = o.Orchestrator.from_execution_parameters(params)
    return orchestrator.execute(**params.operation_parameters)


class GenerateArgParser:
    def __init__(self):
        parser = ArgumentParser()
        subparsers = parser.add_subparsers(dest='operation')
        parser.add_argument('-p', '--parallelism', choices=p.get_parallelism_strings(), nargs='?',
                            const=p.LEVEL_0, default=p.LEVEL_0,
                            help='Selects which level of parallelism to execute the molecule generation with')
        parser.add_argument('-c', '--chunk', '--chunk_size', type=int, default=None, dest='chunk_size',
                            help='The number of generated records to buffer before saving them to the repository. Default is specified in config.py')

        length_parser = ArgumentParser(add_help=False)
        length_parser.add_argument('-l', '--length', '--peptide_length', type=int, choices=[3, 4, 5],
                                         required=True, dest='peptide_length',
                                         help='The number of monomers to assemble into a peptide.')

        sidechain_parser = subparsers.add_parser(g.SidechainModifier.STRING)
        sidechain_parser.set_defaults(func=execute)

        monomer_parser = subparsers.add_parser(g.MonomerGenerator.STRING)
        monomer_parser.set_defaults(func=execute)

        inter_reaction_parser = subparsers.add_parser(g.InterMolecularReactionGenerator.STRING)
        inter_reaction_parser.set_defaults(func=execute)

        intra_reaction_parser = subparsers.add_parser(g.IntraMolecularReactionGenerator.STRING)
        intra_reaction_parser.set_defaults(func=execute)

        peptide_plan_parser = subparsers.add_parser(g.PeptidePlanGenerator.STRING, parents=[length_parser])
        peptide_plan_parser.add_argument('-n', '--num', '--num_peptides', type=int, required=True, dest='num_peptides',
                                         help='The number of peptides to generate.')
        peptide_plan_parser.set_defaults(func=execute)

        peptide_parser = subparsers.add_parser(g.PeptideGenerator.STRING, parents=[length_parser])
        peptide_parser.set_defaults(func=execute)

        template_peptide_parser = subparsers.add_parser(g.TemplatePeptideGenerator.STRING, parents=[length_parser])
        template_peptide_parser.set_defaults(func=execute)

        macrocycle_parser = subparsers.add_parser(g.MacrocycleGenerator.STRING, parents=[length_parser])
        macrocycle_parser.set_defaults(func=execute)

        conformer_parser = subparsers.add_parser(g.ConformerGenerator.STRING, parents=[length_parser])
        conformer_parser.set_defaults(func=execute)

        args = parser.parse_args()
        self.return_val = args.func(args)


if __name__ == "__main__":
    generate_parser = GenerateArgParser()
