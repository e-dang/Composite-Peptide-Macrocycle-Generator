"""
Written by Eric Dang.
github: https://github.com/e-dang
email: edang830@gmail.com
"""

import argparse
import os
import sys
from collections import OrderedDict
from logging import INFO

from rdkit import Chem

from macrocycles.config import (CHAIN_MAP_NUM, METHYL, REGIO_DOC_TYPE,
                                REGIO_INPUT_COL, REGIO_INPUT_DIR,
                                REGIO_OUTPUT_COL, REGIO_OUTPUT_DIR,
                                RXN_MAP_NUM)
from macrocycles.exceptions import AtomSearchError, InvalidSmilesString
from macrocycles.utils import (Base, Flags, IOPaths, MongoParams,
                               create_logger, set_flags, test_valid_smiles)

LOGGER = create_logger(__name__, INFO)


class RegioIsomerEnumerator(Base):
    """
    Class for enumerating the regioisomers of side chains undergoing EAS reactions. Does not consider how likely the
    regioisomer is to occur, rather it only enumerates all atom that satisfy the specified conditions. Inherits from
    Base.

    Attributes:
        side_chains (list): Contains the different side chains for generating regioisomers.
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
        f_in = [os.path.join(REGIO_INPUT_DIR, file) for file in kwargs['f_in']] if 'f_in' in kwargs else ['']
        f_out = os.path.join(REGIO_OUTPUT_DIR, kwargs['f_out']) if 'f_out' in kwargs else ''
        mongo_params = kwargs['mongo_params'] if 'mongo_params' in kwargs else MongoParams(
            REGIO_INPUT_COL, REGIO_DOC_TYPE, REGIO_OUTPUT_COL)
        mongo_params = mongo_params if not no_db else None
        super().__init__(IOPaths(f_in, f_out), mongo_params, LOGGER, input_flags, output_flags)

        # data
        self.side_chains = []

    def load_data(self):
        """
        Overloaded method for loading input data.

        Returns:
            bool: True if successful.
        """

        try:
            self.side_chains = super().load_data()[0]  # should be a single item list, access only item
        except IndexError:
            self.logger.exception('Check MongoParams contains correct number of input_cols and input_types. '
                                  f'self.mongo_params = {self.mongo_params}')
        except Exception:
            self.logger.exception('Unexpected exception occured.')
        else:
            return True

        return False

    def enumerate_regioisomers(self):
        """
        Main driver of class functionality. Calls self.set_atom_map_nums() on each molecule in self.side_chains and
        passes the results to self.accumulate_mols() for proper formatting.

        Returns:
            bool: True if successful.
        """

        try:
            for side_chain in self.side_chains:
                if METHYL in side_chain['modifications']:
                    regioisomers = self.set_atom_map_nums(Chem.MolFromSmiles(side_chain['smiles']))
                    self.accumulate_mols(regioisomers, side_chain)
        except AtomSearchError:
            self.logger.exception()
        except Exception:
            self.logger.exception('Unexpected exception occured.')
        else:
            self.logger.info(f'Successfully generated {len(self.result_data)} regioisomers')
            return True

        return False

    def set_atom_map_nums(self, side_chain):
        """
        Enumerate all possible regioisomer locations on each sidechain for EAS reactions and all possible attachment
        points of side chain to peptide backbone. Possible regioisomer locations are determined by:
            Carbon - has at least 1 hydrogen, is aromatic, and atom is not designated as peptide backbone attachment
                point
            Nitrogen, Oxygen, Sulfer - has at least 1 hydrogen

        Args:
            side_chain (rdkit Mol): The rdkit Mol representation of the side chain

        Returns:
            dict: Keys are SMILES strings, and values are reacting atom index. Each SMILES string is a different
                regioisomer or has a different peptide backbone attachment point.
        """

        regioisomers = {}

        # attachment point at terminal end of alkyl chain
        matches = side_chain.GetSubstructMatches(Chem.MolFromSmarts('[CH3]*'), useChirality=False)
        for pairs in matches:

            # assign atom map number for chain attachment carbon
            for atom_idx in pairs:
                atom = Chem.Mol.GetAtomWithIdx(side_chain, atom_idx)
                if atom.GetSymbol() == 'C' and Chem.Atom.GetTotalNumHs(atom) == 3:
                    atom.SetAtomMapNum(CHAIN_MAP_NUM)
                    old_atom_idx = atom_idx
                    break
            else:
                raise AtomSearchError(f'Couldn\'t find attchment point of: {Chem.MolToSmiles(side_chain)}')

            # set atom map numbers for every eligible EAS reacting atom
            for atom in side_chain.GetAtoms():

                # determine atom eligibility
                if (atom.GetSymbol() == 'C' and atom.GetAtomMapNum() != CHAIN_MAP_NUM and atom.GetTotalNumHs() != 0
                        and atom.GetIsAromatic()) or (atom.GetSymbol() in ['N', 'O', 'S'] and atom.GetTotalNumHs() > 0):
                    atom.SetAtomMapNum(RXN_MAP_NUM)
                else:
                    continue

                # format SMILES string
                atom_idx = atom.GetIdx()
                mol_str = Chem.MolToSmiles(side_chain)
                ind = mol_str.find(f'[CH3:{CHAIN_MAP_NUM}]') + 7  # length of atom map string
                mol_str = mol_str[:ind] + f'[*:{CHAIN_MAP_NUM}]' + mol_str[ind:]
                mol_str = mol_str.replace(f'[CH3:{CHAIN_MAP_NUM}]', 'C')

                atom.SetAtomMapNum(0)   # reset reacting atom map number

                try:
                    test_valid_smiles(mol_str)
                except InvalidSmilesString:
                    self.logger.exception(
                        f'Invalid SMILES string created for side chain {Chem.MolFromSmiles(side_chain)}')
                else:
                    regioisomers[mol_str] = atom_idx

            side_chain.GetAtomWithIdx(old_atom_idx).SetAtomMapNum(0)  # reset attachment atom map number

        return regioisomers

    def accumulate_mols(self, regioisomers, side_chain):
        """
        Stores all data associated with the different regioisomers into a dictionary and appends it so self.result_data.

        Args:
            regioisomers (dict): A dictionary containing the regioisomer's atom mapped SMILES strings as keys and the
                reacting atom index as the corresponding value.
            side_chain (dict): The associated data of the side chain that the regioisomers are derived from.
        """

        for i, (smarts, rxn_atom_idx) in enumerate(regioisomers.items()):
            doc = OrderedDict([('ID', 'r' + str(i) + side_chain['ID']),
                               ('type', 'side_chain'),
                               ('smarts', smarts),
                               ('side_chain', side_chain),
                               ('rxn_atom_idx', rxn_atom_idx)])
            self.result_data.append(doc)


def main():
    """
    Driver function. Parses arguments, constructs class, and performs operations on data.
    """

    parser = argparse.ArgumentParser(description='Enumerates all regioisomers for each side chain using atom map '
                                     'numbers. The atom map numbers associated with the reacting atom and the atom '
                                     'connecting the side chain to the rest of the peptide is defined in config.py.')
    parser.add_argument('input', choices=['json', 'txt', 'mongo', 'sql'],
                        help='Specifies the format that the input data is in.')
    parser.add_argument('output', nargs='+', choices=['json', 'txt', 'mongo', 'sql'],
                        help='Specifies what format to output the result data in.')
    parser.add_argument('--f_in', dest='f_in', required='json' in sys.argv or 'txt' in sys.argv,
                        help='The input file relative to default input directory defined in config.py.')
    parser.add_argument('--f_out', dest='f_out', default='regioisomers',
                        help='The output file relative to the default output directory defined in config.py.')
    parser.add_argument('--no_db', dest='no_db', action='store_true',
                        help='Turns off default connection that is made to the database.')

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

    # create class and perform operations
    enumerator = RegioIsomerEnumerator(f_in=f_in, f_out=args.f_out, input_flags=input_flags,
                                       output_flags=output_flags, no_db=args.no_db)
    if enumerator.load_data() and enumerator.enumerate_regioisomers():
        return enumerator.save_data()

    return False


if __name__ == '__main__':
    main()
