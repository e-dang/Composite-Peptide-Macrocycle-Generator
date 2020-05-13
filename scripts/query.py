from argparse import ArgumentParser

from rdkit import Chem

import cpmg.operations as ops
import cpmg.repository as r
from cpmg.ranges import Key, WholeRange


FIND_OPERATION = 'find'


def execute(args):
    if args.operation == FIND_OPERATION:
        finder = ops.Finder()
        return finder.execute(args)


def print_database(args):
    print(r.create_repository_from_string(args.repo))

    return True


def remove_dataset(args):
    repo = r.create_repository_from_string(args.repo)
    for dataset in args.dataset:
        full_path = '/'.join([args.repo, dataset])
        print(f'Removing dataset {full_path}...')
        repo.remove_group(dataset)


def execute_deactivate(args):
    repo = r.create_repository_from_string(args.repo)
    repo_mols = {mol.kekule: mol._id for mol in repo.load()}

    file_mols = []
    reader = Chem.SDMolSupplier(args.file)
    for mol in reader:
        mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
        Chem.Kekulize(mol)
        file_mols.append(Chem.MolToSmiles(mol, kekuleSmiles=True))

    if args.invert:
        deactivated_mols = set(repo_mols.keys()).difference(file_mols)
    else:
        deactivated_mols = file_mols

    deactivated_ids = [repo_mols[kekule] for kekule in deactivated_mols]
    return repo.deactivate_records(Key(deactivated_ids))


def execute_activate(args):
    repo = r.create_repository_from_string(args.repo)
    if args.all:
        ids = [record._id for record in repo.load_inactive_records(Key(WholeRange()))]
        return repo.activate_records(Key(ids))


class QueryArgParser:
    def __init__(self):
        parser = ArgumentParser()
        subparsers = parser.add_subparsers(dest='operation')

        print_parser = subparsers.add_parser('print')
        print_parser.add_argument('repo', choices=r.get_all_repository_strings(), nargs='?',
                                  const=r.CPMGRepository.STRING,
                                  default=r.CPMGRepository.STRING,
                                  help='Prints the number of records and structure if applicable in the repository for the specified type. Default type is all.')
        print_parser.set_defaults(func=print_database)

        find_parser = subparsers.add_parser(FIND_OPERATION)
        find_parser.add_argument('repo', choices=r.get_all_repository_strings(),
                                 help='The repository to search.')
        find_parser.add_argument('-a', '--all', action='store_true',
                                 help='Flag that causes all records in the specified repo to be loaded.')
        find_parser.add_argument('-k', '--kekule', type=str, nargs='+',
                                 help='Load and print the records corresponding to the specified kekule strings.')
        find_parser.add_argument('-f', '--filter', type=str,
                                 help='The filter criteria as a python dictionary, where the key is the field on the record and the value is the value of that field in desired record(s). File must be a .sdf file')
        find_parser.add_argument('-p', '--projection', type=str, nargs='+',
                                 help='Only print the specified field from the retrieved records.')
        find_parser.add_argument('-i', '--inactives', action='store_true',
                                 help='Flag that selects inactive records rather than active records in the given repository that meet the search criteria.')
        find_parser.add_argument('-c', '--chunk', '--chunk_size', type=int, const=10, default=10, nargs='?', dest='chunk_size',
                                 help='The number of records to display on screen at a single time.')
        find_parser.set_defaults(func=execute)

        remove_parser = subparsers.add_parser('remove')
        remove_parser.add_argument('repo', type=str,
                                   help='The repository that contains the dataset to be removed.')
        remove_parser.add_argument('-d', '--dataset', type=str, nargs='+', help='The dataset to remove.')
        remove_parser.set_defaults(func=remove_dataset)

        deactivate_parser = subparsers.add_parser('deactivate')
        deactivate_parser.add_argument('repo', type=str,
                                       help='The repository that contains the records to be deactivated.')
        deactivate_parser.add_argument('-f', '--file', type=str, required=True,
                                       help='The file containing the structures to inactivate.')
        deactivate_parser.add_argument('-i', '--invert', action='store_true',
                                       help='Causes the structures in the file to remain active while all records not matched to structures in the file are deactivated.')
        deactivate_parser.set_defaults(func=execute_deactivate)

        activate_parser = subparsers.add_parser('activate')
        activate_parser.add_argument('repo', type=str,
                                     help='The repository that contains the records to be deactivated.')
        activate_parser.add_argument('-a', '--all', action='store_true',
                                     help='Switch that causes all records within the repository to be activated.')
        activate_parser.set_defaults(func=execute_activate)

        args = parser.parse_args()
        self.return_val = args.func(args)


if __name__ == "__main__":
    query_parser = QueryArgParser()
