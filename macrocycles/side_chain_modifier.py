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

from macrocycles.config import (COL1, CONN_MAP_NUM, HETERO_MAP_NUM,
                                SCM_DOC_TYPE, SCM_INPUT_DIR, SCM_OUTPUT_DIR,
                                Connections)
from macrocycles.exceptions import MissingMapNumberError
from macrocycles.utils import (Base, IOPaths, Mongo, Smiles,
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

    def __init__(self, f_in, f_out, input_flags, output_flags):
        """
        Constructor.

        Args:
            f_in (str): The file name containing the input data, if input data has been specified to be in file format.
            f_out (str): The file name that the result_data will be written to, if specified to do so.
            input_flags (Flags): A namedtuple containing the flags indicating which format to get input data from.
            output_flags (Flags): A namedtuple containing the flags indicating which format to output data to.
        """

        # I/O
        f_in = os.path.join(SCM_INPUT_DIR, f_in)
        f_out = os.path.join(SCM_OUTPUT_DIR, f_out)
        super().__init__(IOPaths(f_in, f_out), mongo_setup=Mongo(COL1, SCM_DOC_TYPE),
                         logger=LOGGER, input_flags=input_flags, output_flags=output_flags)

        # data
        self.parent_side_chains = []
        self.connections = [Connections(con, mod) for con, mod in [(f'[CH3:{CONN_MAP_NUM}]', [0, 3]),
                                                                   (f'[CH3][CH2:{CONN_MAP_NUM}]', [1, 3]),
                                                                   (f'[CH3][CH2][CH2:{CONN_MAP_NUM}]', [2, 3])]]

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
            return True

        return False

    def load_data(self):
        """
        Overloaded method for loading input data.

        Returns:
            bool: True if successful.
        """

        try:
            self.parent_side_chains = super().load_data()
        except Exception:
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
            set: A set of unique SMILES strings representing the side chain with alkyl chains attached at different
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
            mols (iterable): A set containing the unique SMILES strings of the modified side chain.
            parent (rdkit Mol): The parent side chain from which the modified side chain was derived.
            modifications (iterable): A list containing ids of the modifications that were made to the side chain.
            group (str): The group name that the parent side chain came from.
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
                                     'chains to all elgible positions on the side chain. Alkyl chains include '
                                     "methyl, ethyl, and propyl. The last word in the input files' name dictates which "
                                     'group the side chain belongs to.')
    parser.add_argument('input', choices=['json', 'txt', 'mongo', 'sql'],
                        help='Specifies the format that the input data is in.')
    parser.add_argument('output', nargs='+', choices=['json', 'txt', 'mongo', 'sql'],
                        help='Specifies what format to output the result data in.')
    parser.add_argument('--f_in', dest='f_in', required='json' in sys.argv or 'txt' in sys.argv,
                        help='The input file relative to default input directory defined in config.py.')
    parser.add_argument('--f_out', dest='f_out', default='side_chains',
                        help='The output file relative to the default output directory defined in config.py.')

    args = parser.parse_args()

    # check for proper file specifications
    if args.input in ['json', 'txt']:
        extension = args.f_in.split('.')[-1]
        if args.input != extension:
            LOGGER.error('File extension of the input file does not match the specified format')
            raise OSError('File extension of the input file does not match the specified format')

    # configure I/O
    input_flags, output_flags = set_flags(args.input, args.output)
    f_in = args.f_in if args.input in ['json', 'txt'] else '/'

    # create class and perform operations
    modifier = SideChainModifier(f_in, args.f_out, input_flags=input_flags, output_flags=output_flags)
    if modifier.load_data() and modifier.diversify():
        return modifier.save_data()

    return False


if __name__ == '__main__':
    main()
