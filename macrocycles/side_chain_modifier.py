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

from macrocycles.config import (CONN_MAP_NUM, HETERO_MAP_NUM,
                                SCM_DOC_TYPE, SCM_INPUT_DIR, SCM_OUTPUT_DIR, SCM_INPUT_COL, SCM_OUTPUT_COL, CONNECTIONS)
from macrocycles.exceptions import MissingMapNumberError
from macrocycles.utils import (Base, IOPaths, Flags, MongoParams, Smiles,
                               create_logger, set_flags)

LOGGER = create_logger(name=__name__, level=INFO)


class SideChainModifier(Base):
    """
    Class for attaching varying length alkyl chains / attachment points to side chain. Inherits from Base.

    Attributes:
        parent_side_chains (list): Contains the parent_side_chains and associated data as dictionaries.
        connections (list): Contains the atom mapped SMARTS strings of the alkyl attachment chain and the corresponding
            modification array as a Connections named tuple.
    """

    def __init__(self, logger=LOGGER, input_flags=Flags(False, False, True, False),
                 output_flags=Flags(False, False, False, False), no_db=False, **kwargs):
        """
        Constructor.

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
        f_in = [os.path.join(SCM_INPUT_DIR, file) for file in kwargs['f_in']] if 'f_in' in kwargs else ['']
        f_out = os.path.join(SCM_OUTPUT_DIR, kwargs['f_out']) if 'f_out' in kwargs else ''
        mongo_params = kwargs['mongo_params'] if 'mongo_params' in kwargs else MongoParams(
            SCM_INPUT_COL, SCM_DOC_TYPE, SCM_OUTPUT_COL)
        mongo_params = mongo_params if not no_db else None
        super().__init__(IOPaths(f_in, f_out), mongo_params, LOGGER, input_flags, output_flags)

        # data
        self.parent_side_chains = []
        self.connections = CONNECTIONS

    def diversify(self):
        """
        Main driver of class functionality. Calls self.alternate_connection_point() on each
        molecule in self.parent_side_chains with each connection type in self.connections and calls
        self.accumulate_mols() on the resulting data.

        Returns:
            bool: True if successful.
        """

        for doc in self.parent_side_chains:
            unique_mols = []
            for connection, modification in self.connections:
                try:
                    unique_mols.extend(self.alternate_connection_point(doc['smiles'], connection))
                except MissingMapNumberError:
                    self.logger.exception(f'Connection missing map numbers! connection: {connection}')
                    break
            else:
                self.accumulate_mols(unique_mols, doc, modification)
                continue
            break
        else:
            self.logger.info(f'Successfully modified parent side chains into {len(self.result_data)} side chains')
            return True

        return False

    def load_data(self):
        """
        Overloaded method for loading input data.

        Returns:
            bool: True if successful.
        """

        try:
            self.parent_side_chains = super().load_data()[0]  # single item list is returned, access only item
        except IndexError:
            self.logger.exception('Check MongoParams contains correct number of input_cols and input_types. '
                                  f'self.mongo_params = {self.mongo_params}')
            return False
        except Exception:
            self.logger.exception('Unexcepted exception occured.')
            return False
        else:
            return True

    def alternate_connection_point(self, mol, connection):
        """
        Creates a set of new molecules by attaching an alkyl chain (which becomes the attachment point to peptide
        backbone) to every eligble position on the side chain. Eligiblity of an atom is defined as:
            Carbon - Must have 1 or 2 hydrogens
            Nitrogen, Oxygen, Sulfur - Must have > 0 hydrogens

        Args:
            mol (rdkit Mol / str): The side chain molecule.
            connection (rdkit Mol / str): The atom mapped alkyl attachment chain. If in string format, must be a valid
                SMARTS string.

        Raises:
            MissingMapNumberError: If connection is not atom mapped.

        Returns:
            list: A list of unique SMILES strings representing the side chain with alkyl chains attached at different
                positions.
        """

        mols = set()
        attached = set()

        # convert strings to rdkit Mol
        if isinstance(connection, str):
            connection = Chem.MolFromSmarts(connection)
        if isinstance(mol, str):
            mol = Chem.MolFromSmiles(mol)

        # check if connecting atom is atom mapped
        map_nums = [atom.GetAtomMapNum() for atom in connection.GetAtoms()]
        if CONN_MAP_NUM not in map_nums:
            raise MissingMapNumberError('Need to specifiy connecting atom with atom map number')

        # make attachment at each atom
        for atom in mol.GetAtoms():

            # detetmine atom eligibility
            atom_idx = None
            found = False
            if atom.GetSymbol() == 'C' and atom.GetTotalNumHs() != 0 and atom.GetTotalNumHs() < 3:
                atom.SetAtomMapNum(HETERO_MAP_NUM)
                atom_idx = atom.GetIdx()
                attached.add(atom_idx)
                found = True
            elif atom.GetSymbol() in ('N', 'O', 'S') and atom.GetTotalNumHs() != 0:
                atom.SetAtomMapNum(HETERO_MAP_NUM)
                atom_idx = atom.GetIdx()
                attached.add(atom_idx)
                found = True
            elif atom.GetIdx() in attached or not found:
                continue

            # prepare for attachment
            combo = Chem.RWMol(Chem.CombineMols(mol, connection))

            # find reacting atoms on combo mol and reset atom map numbers
            mol_atom = None
            conn_atom = None
            for combo_atom in combo.GetAtoms():
                if combo_atom.GetAtomMapNum() == HETERO_MAP_NUM:
                    mol_atom = combo_atom.GetIdx()
                    combo_atom.SetAtomMapNum(0)
                    Chem.Mol.GetAtomWithIdx(mol, atom_idx).SetAtomMapNum(0)
                elif combo_atom.GetAtomMapNum() == CONN_MAP_NUM:
                    conn_atom = combo_atom.GetIdx()
                    combo_atom.SetAtomMapNum(0)

            combo.AddBond(mol_atom, conn_atom, order=Chem.rdchem.BondType.SINGLE)

            # fix hydrogen counts
            atom_react = combo.GetAtomWithIdx(mol_atom)
            if atom_react.GetSymbol() in ('N', 'O', 'S'):
                atom_react.SetNumExplicitHs(0)
            elif atom_react.GetSymbol() == 'C' and Chem.Atom.GetTotalNumHs(atom_react) > 0:
                atom_react.SetNumExplicitHs(Chem.Atom.GetTotalNumHs(atom_react) - 1)

            try:
                Chem.SanitizeMol(combo)
                Chem.Kekulize(combo)
            except ValueError:
                self.logger.exception(f'Sanitize error! Side Chain: {Chem.MolToSmiles(mol)}')
            else:
                mols.add(Smiles(Chem.MolToSmiles(combo), Chem.MolToSmiles(combo, kekuleSmiles=True)))
                self.logger.debug(f'Success! Side Chain: {Chem.MolToSmiles(mol)}')

        return list(mols)

    def accumulate_mols(self, unique_mols, parent, modifications):
        """
        Stores all data associated with the modified side chain in a dictionary and appends it to self.result_data.

        Args:
            unique_mols (iterable): A list containing the unique SMILES strings of the modified side chain.
            parent (dict): The assocaited data of parent side chain from which the modified side chain was derived.
            modifications (list): A list containing ids of the modifications that were made to the parent side chain.
        """

        for i, smiles in enumerate(unique_mols):
            doc = OrderedDict([('ID', parent['ID'] + str(i)),
                               ('type', 'side_chain'),
                               ('smiles', smiles.smiles),
                               ('kekule', smiles.kekule),
                               ('modifications', modifications),
                               ('parent', parent)])
            self.result_data.append(doc)


def main():
    """
    Driver function. Parses arguments, constructs class, and performs operations on data.
    """

    parser = argparse.ArgumentParser(description='Creates a unique set of molecules by attaching varying length alkyl '
                                     'chains to all elgible positions on the parent side chain. Alkyl chains include '
                                     'methyl, ethyl, and propyl.')
    parser.add_argument('input', choices=['json', 'txt', 'mongo', 'sql'],
                        help='Specifies the format that the input data is in.')
    parser.add_argument('output', nargs='+', choices=['json', 'txt', 'mongo', 'sql'],
                        help='Specifies what format to output the result data in.')
    parser.add_argument('--f_in', dest='f_in', required='json' in sys.argv or 'txt' in sys.argv,
                        help='The input file relative to default input directory defined in config.py.')
    parser.add_argument('--f_out', dest='f_out', default='side_chains',
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
    modifier = SideChainModifier(f_in=f_in, f_out=args.f_out, input_flags=input_flags,
                                 output_flags=output_flags, no_db=args.no_db)
    if modifier.load_data() and modifier.diversify():
        return modifier.save_data()

    return False


if __name__ == '__main__':
    main()
