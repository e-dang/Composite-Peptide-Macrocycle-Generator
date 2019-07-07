"""
Written by Eric Dang.
github: https://github.com/e-dang
email: edang830@gmail.com
"""

import argparse
import multiprocessing
import os
import sys
from collections import OrderedDict
from copy import copy
from itertools import islice, product
from logging import INFO
from time import time

from rdkit import Chem

from macrocycles.config import (CAPACITY, CARBON_MAP_NUM, NITROGEN_MAP_NUM,
                                PG_DOC_TYPE, PG_INPUT_COL, PG_INPUT_DIR,
                                PG_OUTPUT_COL, PG_OUTPUT_DIR)
from macrocycles.utils import (Base, Flags, IOPaths, MongoParams,
                               create_logger, ranges, set_flags)

LOGGER = create_logger(__name__, INFO)


class PeptideGenerator(Base):
    """
    Class for joining monomers together to form a peptide. Inherits from Base.

    Attributes:
        monomers (list): Contains the monomers and the associated data.
    """

    def __init__(self, logger=LOGGER, input_flags=Flags(False, False, True, False),
                 output_flags=Flags(False, False, False, False), no_db=False, **kwargs):
        """
        Initializer.

        Args:
            input_flags (Flags): A namedtuple containing the flags indicating which format to get input data from.
            output_flags (Flags): A namedtuple containing the flags indicating which format to output data to.
            no_db (bool): If True, ensures that no default connection is made to the database. Defaults to False.
            **kwargs:
                f_in (str): The file name(s) containing the input data, if input data has been specified to be retrieved
                    from a file.
                f_out (str): The file name that the result_data will be written to, if specified to do so.
                mongo_params (MongoParams): A namedtuple containing the collection name(s) and value(s) held within the
                    'type' field of the documents to be retrieved, as well as the output collection name.
                no_db (bool): If True, ensures that no default connection is made to the database.
        """

        # I/O
        f_in = [os.path.join(PG_INPUT_DIR, file) for file in kwargs['f_in']] if 'f_in' in kwargs else ['']
        f_out = os.path.join(PG_OUTPUT_DIR, kwargs['f_out']) if 'f_out' in kwargs else ''
        mongo_params = kwargs['mongo_params'] if 'mongo_params' in kwargs else MongoParams(
            PG_INPUT_COL, PG_DOC_TYPE, PG_OUTPUT_COL)
        mongo_params = mongo_params if not no_db else None
        super().__init__(IOPaths(f_in, f_out), mongo_params, LOGGER, input_flags, output_flags)

        # data
        self.monomers = []

    def load_data(self):
        """
        Overloaded method for loading input data.

        Returns:
            bool: True if successful.
        """

        try:
            self.monomers = list(super().load_data()[0])
        except ValueError:
            self.logger.exception('Check MongoParams contains correct number of input_cols and input_types. '
                                  f'self.mongo_params = {self.mongo_params}')
        except Exception:
            self.logger.exception('Unexcepted exception occured.')
        else:
            return True

        return False

    def generate_peptides(self, length, num_peptides=None, num_jobs=1, job_id=1):
        """
        Top level function that generates the cartesian product of monomers in self.monomers. create_peptide() is then
        called in a parallel fashion on each tuple of monomers is the cartesian product. The resulting peptide and
        monomers are then passed to self.accumulate_data(). If self.result_data is about to exceed the capacity of a
        python list, an automatic call to self.save_data() is made and self.result_data is cleared.

        Args:
            length (int): The number of monomers per peptide.
            num_peptides (int, optional): The number of peptides to make. If None, then all peptides in the cartesian
                product as produced. Defaults to None.
            num_jobs (int, optional): The number of jobs submitted in a job array. Used for running this method on a
                super computer. Does not do anything if num_peptides is specified. Defaults to 1.
            job_id (int, optional): The job ID of the specific job in the job array running this function. Used for
                running this method on a super computer. Does not do anything if num_peptides is specified. Defaults
                to 1.

        Returns:
            bool: True if successful.
        """

        def save_data(self, num_peptides, job_id):
            """
            Overloaded variant for saving data with different file names that depend on how many times this method has
            been called and the job_id if specified. Necessary since python lists may not be able to hold the number of
            peptides that are generated.

            Args:
                num_peptides (int): Used to determine whether num_peptides was specified in call to generate_peptides().
                job_id (int): The job_id passed to generate_peptides().
            """

            # intialize counting variable
            if 'count' not in self.__dict__.keys():
                self.count = 0

            # format file name based on job_id and how many times this function has been called
            old_fp = copy(self.fp_out)
            self.fp_out = self.fp_out.split('.')
            if not num_peptides:
                self.fp_out[0] += '_' + str(job_id)
            self.fp_out[0] += '_' + str(self.count)
            self.fp_out = '.'.join(self.fp_out)
            self.count += 1

            self.save_data()
            self.fp_out = old_fp

        try:

            # determine how to run method
            if num_peptides:
                monomers = islice(product(self.monomers, repeat=length), num_peptides)
            else:
                start, stop = ranges(len(self.monomers) ** length, num_jobs)[job_id - 1]
                monomers = islice(product(self.monomers, repeat=length), start, stop)

            # perform peptide generation in parallel
            with multiprocessing.Pool() as pool:
                result = pool.map_async(PeptideGenerator.create_peptide, monomers, chunksize=25000)
                for peptide, monomer in result.get():
                    self.accumulate_data(peptide, monomer)

                    if len(self.result_data) > CAPACITY:
                        save_data(self, num_peptides, job_id)
                        self.result_data = []
                # for peptide, monomer in pool.imap_unordered(self.connect_monomers, monomers, chunksize=25000):
                # if peptide is not None or monomer is not None:
                # pass
                # self.accumulate_data(peptide, monomer)
        except Exception:
            self.logger.exception('Unexpected exception occured.')
        else:
            return True

        return False

    @staticmethod
    def create_peptide(monomers, logger=None):
        """
        Static class method. Connects monomers together creating a peptide and verifies that the peptide is composed of
        at least one required monomer.

        Args:
            monomers (list): Contains the monomer documents.
            logger (Logger, optional): A logger used to log errors in the method. Defaults to None.

        Returns:
            dict: A dictionary containing the peptide SMILES string as a key and a list containing the kekule SMILES
                string and N-terminal monomer's backbone type as a value.
            list: The list of monomer documents used to create the peptide.
        """

        def is_valid(monomers):
            """
            Determines the validity of the peptide.

            Args:
                monomers (list): Contains the monomer documents.

            Returns:
                bool: True if at least one monomer is required.
            """

            for monomer in monomers:
                if monomer['required']:
                    return True

            return False

        # check taht peptide is valid before constructing
        if not is_valid(monomers):
            return None, None

        pattern_dict = {'alpha': Chem.MolFromSmarts('NCC(=O)O'),
                        'beta2': Chem.MolFromSmarts('NCCC(=O)O'),
                        'beta3': Chem.MolFromSmarts('NCCC(=O)O')}

        # begin conneting each monomer in monomers
        for i, monomer_doc in enumerate(monomers):

            monomer = Chem.MolFromSmiles(monomer_doc['smiles'])

            # start peptide with first monomer
            if i == 0:
                peptide = monomer
                previous_type = monomer_doc['backbone']
                continue

            peptide_match = peptide.GetSubstructMatches(pattern_dict[previous_type])
            monomer_match = monomer.GetSubstructMatches(pattern_dict[monomer_doc['backbone']])
            previous_type = monomer_doc['backbone']

            # set atom map number for nitrogen connection point of monomer
            for pairs in monomer_match:
                for atom_idx in pairs:
                    atom = monomer.GetAtomWithIdx(atom_idx)
                    if atom.GetSymbol() == 'N' and atom.GetTotalNumHs() != 0:
                        atom.SetAtomMapNum(NITROGEN_MAP_NUM)
                        monomer_old_attach = atom
                        break

            # set atom map number for carboxyl carbon connection point of peptide and remove carboxyl oxygen
            peptide = Chem.RWMol(peptide)
            for pairs in peptide_match:
                for atom_idx in pairs:
                    atom = peptide.GetAtomWithIdx(atom_idx)
                    if atom.GetSymbol() == 'O' and atom.GetTotalNumHs() == 1:
                        carboxyl_atom = atom_idx
                    elif atom.GetSymbol() == 'C' and atom.GetTotalNumHs() == 0 and \
                            atom.GetHybridization() == Chem.rdchem.HybridizationType.SP2:
                        atom.SetAtomMapNum(CARBON_MAP_NUM)
                        pep_old_attach = atom

            peptide.RemoveAtom(carboxyl_atom)

            # connect peptide and monomer
            try:
                peptide = Base.merge(peptide, monomer, CARBON_MAP_NUM, NITROGEN_MAP_NUM)
            except ValueError:
                if logger:
                    logger.exception(f'Sanitize Error! monomer = {Chem.MolToSmiles(monomer)}, '
                                     f'peptide = {Chem.MolToSmiles(peptide)}')
            else:
                pep_old_attach.SetAtomMapNum(0)
                monomer_old_attach.SetAtomMapNum(0)

        smiles = Chem.MolToSmiles(peptide)
        Chem.Kekulize(peptide)
        result = {smiles: [Chem.MolToSmiles(peptide, kekuleSmiles=True), monomer_doc['backbone']]}
        return result, monomers

    def accumulate_data(self, peptide, monomers):
        """
        Stores all data associated with the peptide in a dictionary and appends it to self.result_data.

        Args:
            peptide (dict): The dictionary containing the peptide's SMILES string as keys and a list containing the
                kekule SMILES and the N-terminal monomer's backbone type as values.
            monomers (dict): The dictionary containing the associated monomer data.
        """

        monomer_ids = [monomer['ID'] for monomer in monomers]
        pep_id = ''.join(monomer_ids)
        for smiles, vals in peptide.items():
            doc = OrderedDict([('ID', pep_id),
                               ('type', 'peptide'),
                               ('smiles', smiles),
                               ('kekule', vals[0]),
                               ('length', len(monomers)),
                               ('monomers', monomer_ids),
                               ('N_term', vals[1])])
        self.result_data.append(doc)


def main():
    """
    Driver function. Parses arguments, constructs class, and performs operations on data.
    """

    parser = argparse.ArgumentParser(
        description='Connects specified number of monomers to form a peptide. Takes input file(s) containing '
        'monomer SMILES strings and outputs to a json file the peptides as SMILES strings. Input and output are json '
        'files because this script will run on a super computer')
    parser.add_argument('length', type=int, help='The length of the peptide in monomers.')
    parser.add_argument('input', choices=['json', 'txt', 'mongo', 'sql'],
                        help='Specifies the format that the input data is in.')
    parser.add_argument('output', nargs='+', choices=['json', 'txt', 'mongo', 'sql'],
                        help='Specifies what format to output the result data in.')
    parser.add_argument('--f_in', dest='f_in', required='json' in sys.argv or 'txt' in sys.argv,
                        help='The input file containing backbone data relative to default input directory defined '
                        'in config.py.')
    parser.add_argument('--f_out', dest='f_out', default='custom',
                        help='The output file relative to the default output directory defined in config.py.')
    parser.add_argument('--no_db', dest='no_db', action='store_true',
                        help='Turns off default connection that is made to the database.')
    parser.add_argument('--num_peptides', type=int, dest='num_peptides', help='The number of peptides to generate')
    parser.add_argument('--num_jobs', dest='num_jobs', default=1, type=int,
                        help='The number of jobs to run on a job array.')
    parser.add_argument('--job_id', dest='job_id', default=1, type=int, help='The job number or ID.')

    args = parser.parse_args()

    # check for proper file specifications
    if args.input in ['json', 'txt']:
        extension = args.f_in.split('.')[-1]
        if args.input != extension:
            LOGGER.error('File extension of the input file does not match the specified format')
            raise OSError('File extension of the input file does not match the specified format')

    # configure I/O
    input_flags, output_flags = set_flags(args.input, args.output)
    f_in = [args.f_in] if args.input in ['json', 'txt'] else ['']

    generator = PeptideGenerator(f_in=f_in, f_out=args.f_out, input_flags=input_flags,
                                 output_flags=output_flags, no_db=args.no_db)
    # t = time()
    if generator.load_data() and generator.generate_peptides(args.length, args.num_peptides,
                                                             args.num_jobs, args.job_id):
        # print(len(generator.result_data))
        # print(generator.result_data[0])
        # print(time() - t)
        return generator.save_data()

    return False


if __name__ == '__main__':
    main()
