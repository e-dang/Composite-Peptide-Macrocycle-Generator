import ast
from functools import partial

from rdkit import Chem
from rdkit.Chem import AllChem

import cpmg.repository as repo
import cpmg.utils as utils
from cpmg.ranges import Key


def conformer_to_pdb(key, filepath):
    conformers = list(repo.create_conformer_repository().load(key))
    for conformer in conformers:
        Chem.MolToPDBFile(conformer.mol, utils.rotate_file(filepath))

    return True


class CalcMW:
    STRING = 'mw'

    def __call__(self, mol):
        return (AllChem.CalcExactMolWt(Chem.Mol(mol.binary)), )


class CalcRB:
    STRING = 'rb'

    def __call__(self, mol):
        return (AllChem.CalcNumRotatableBonds(Chem.Mol(mol.binary)), )


class CalcTPSA:
    STRING = 'tpsa'

    def __call__(self, mol):
        return (AllChem.CalcTPSA(Chem.Mol(mol.binary), includeSandP=True), )


def get_calculation_from_string(string):
    if string == CalcMW.STRING:
        return CalcMW()

    if string == CalcRB.STRING:
        return CalcRB()

    if string == CalcTPSA.STRING:
        return CalcTPSA()


def pca_convex_hull(filepath, output):
    set_to_array = lambda set_data: np.array([np.array(list(point)) for point in set_data])

    points = {}
    with open(filepath, 'r') as file:
        csv_reader = csv.reader(file)
        for i, line in enumerate(csv_reader):
            if i != 0:
                points[tuple(map(float, line[1:3]))] = (line[0], line[-1])

    set_data = set(point for point in points.keys())
    array_data = set_to_array(set_data)

    hull_points = set()

    while len(hull_points) < 10000:
        print(len(hull_points))
        hull = ConvexHull(array_data)
        for simplex in hull.simplices:
            for point in zip(array_data[simplex, 0], array_data[simplex, 1]):
                hull_points.add(tuple(point))

        set_data = set_data - hull_points
        array_data = set_to_array(set_data)

    with open(output, 'w') as file:
        for point in hull_points:
            file.write(points[point][0] + ',' + points[point][1] + '\n')


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
