from argparse import ArgumentParser
import cpmg.repository as r
from cpmg.ranges import Key


def print_database(args):
    print(r.create_repository_from_string(args.repo))

    return True


def execute_find(args):
    if args.kekule is not None:
        return print_records_from_kekule(args.kekule)


def print_records_from_kekule(kekules):
    repo = r.CPMGRepository()
    for record in repo.load(Key(kekules, index='kekule')):
        print(record, '\n')

    return True


def remove_dataset(args):
    repo = r.create_repository_from_string(args.repo)
    for dataset in args.dataset:
        full_path = '/'.join([args.repo, dataset])
        print(f'Removing dataset {full_path}...')
        repo.remove_group(dataset)


class QueryArgParser:
    def __init__(self):
        parser = ArgumentParser()
        subparsers = parser.add_subparsers()

        print_parser = subparsers.add_parser('print')
        print_parser.add_argument('repo', choices=r.get_all_repository_strings(), nargs='?',
                                  const=r.CPMGRepository.STRING,
                                  default=r.CPMGRepository.STRING,
                                  help='Prints the number of records and structure if applicable in the repository for the specified type. Default type is all.')
        print_parser.set_defaults(func=print_database)

        find_parser = subparsers.add_parser('find')
        find_parser.add_argument('-k', '--kekule', type=str, nargs='+',
                                 help='Load and print the records corresponding to the specified kekule strings.')
        find_parser.set_defaults(func=execute_find)

        remove_parser = subparsers.add_parser('remove')
        remove_parser.add_argument('repo', type=str,
                                   help='The repository that contains the dataset to be removed.')
        remove_parser.add_argument('-d', '--dataset', type=str, nargs='+', help='The dataset to remove.')
        remove_parser.set_defaults(func=remove_dataset)

        args = parser.parse_args()
        self.return_val = args.func(args)


if __name__ == "__main__":
    query_parser = QueryArgParser()
