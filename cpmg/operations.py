import ast
from functools import partial

from rdkit import Chem

import cpmg.repository as repo
import macrocycles.utils as utils
from cpmg.ranges import Key


def conformer_to_pdb(key, filepath):
    conformers = list(repo.create_conformer_repository().load(key))
    for conformer in conformers:
        Chem.MolToPDBFile(conformer.mol, utils.file_rotator(filepath))

    return True


class Finder:
    def __init__(self):
        self.records = []

    def execute(self, command_line_args):
        self.repo = repo.create_repository_from_string(command_line_args.repo)
        self.__load_records(command_line_args)
        projection = self.__create_projection_func(command_line_args.projection)
        self.__print_records(projection, command_line_args.chunk_size)

        return True

    def __load_records(self, command_line_args):
        if command_line_args.all:
            self.__load_all_records()
        elif command_line_args.kekule:
            self.__load_records_from_kekule(command_line_args.kekule)
        elif command_line_args.filter:
            self.__load_records_from_filter(ast.literal_eval(command_line_args.filter))

    def __load_all_records(self):
        self.records = self.repo.load()

    def __load_records_from_kekule(self, kekules):
        self.records = self.repo.load(Key(kekules, index='kekule'))

    def __load_records_from_filter(self, filter_criteria):
        def filter_func(x, filter_criteria):
            for key, val in filter_criteria.items():
                try:
                    attr = getattr(x, key)
                except AttributeError:
                    attr = x[key]

                if isinstance(attr, dict) and isinstance(val, dict):
                    return filter_func(attr, val)

                if attr != val:
                    return False

            return True

        if self.records is None:
            self.__load_all_records()

        func = partial(filter_func, filter_criteria=filter_criteria)
        self.records = filter(func, self.records)

    def __create_projection_func(self, projection_criteria):
        if projection_criteria is None:
            return lambda x: x

        def projection_func(record):
            for key in list(record.__dict__):
                if key not in projection_criteria:
                    del record.__dict__[key]

            return record

        return projection_func

    def __print_records(self, projection, chunk_size):
        for i, record in enumerate(map(projection, self.records)):
            print(record, '\n')
            if i % chunk_size == 0 and i != 0:
                user_input = input('')
                if user_input.lower() in ('q', 'quit'):
                    break
