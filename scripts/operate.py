from argparse import ArgumentParser
import cpmg.operations as ops
from cpmg.ranges import Key, WholeRange
import cpmg.models as m

PDB_OPERATION = 'pdb'


def execute(args):
    if args.type == m.Conformer.STRING:
        return execute_conformer_ops(args)


def execute_conformer_ops(args):
    if args.operation == PDB_OPERATION:
        key = generate_key_from_args(args)
        return ops.conformer_to_pdb(key, args.file)


def generate_key_from_args(args):
    if args.all:
        return Key(WholeRange())


class OperateArgParser:
    def __init__(self):
        parser = ArgumentParser()
        subparsers = parser.add_subparsers(dest='type')

        conformer_parser = subparsers.add_parser(m.Conformer.STRING)
        conformer_subparsers = conformer_parser.add_subparsers(dest='operation')
        conformer_pdb_parser = conformer_subparsers.add_parser(PDB_OPERATION)
        conformer_pdb_parser.add_argument('-f', '--file', type=str, required=True,
                                          help='The filepath containing the base filename of where the pdb files will be written.')
        conformer_pdb_parser.add_argument('-a', '--all', action='store_true',
                                          help='Causes all conformers in the repository to be output to pdb files.')
        conformer_pdb_parser.set_defaults(func=execute)

        args = parser.parse_args()
        self.return_val = args.func(args)


if __name__ == "__main__":
    operate_parser = OperateArgParser()
