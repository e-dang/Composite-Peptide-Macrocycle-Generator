"""
Written by Eric Dang.
github: https://github.com/e-dang
email: erickdang@g.ucla.edu
"""

import argparse
from collections import Iterable
from os.path import join

from rdkit import Chem

from base import Base
from database import MolDatabase
from utils import read_mols

# TODO: allow for side chains with methyl group already attached, such as 5-methyl indole.


class SideChainModifier(Base):
    """
    Class for attaching varying length alkyl chains / attachment points to side_chain.

    Attributes:
        _fp_in (list): Contains all filepaths to input files containing the parent side_chain.
        _fp_out (str): The filepath to the output json file.
        _mol_db (MolDatabase): An instance of MolDatabase.
        _groups (list): Contains the groups from which each set of side_chains comes from.
        _connections (list): Contains tuples, each of which contain atom mapped SMARTS strings of the alkyl attachment
            chain and the corresponding modification array.
        _data (list): Contains the dictionaries associated with the modified side_chain.
        _side_chains (list): Contains the parent side_chains as rdkit Mols
    """

    HETERO_MN = 1
    CONN_MN = 2

    def __init__(self, fp_in, fp_out, groups=None, database=None, json_flag=True, db_flag=True):
        """
        Constructor.

        Args:
            fp_in (list): Contains all filepaths to input files containing the parent side_chain.
            fp_out (str): The filepath to the output json file.
            groups (list): Contains the groups from which each set of side_chains comes from.
            database (str, optional): The database to store data in. Defaults to MolDatabase
            json_flag (bool): Determines whether to write data to the json file.
            db_flag (bool): Determines whether to write data to the database.
        """

        # I/O
        super().__init__()
        self.fp_in = [str(self.project_dir / file) for file in fp_in]
        self.fp_out = str(self.project_dir / fp_out)
        self.mol_db = MolDatabase() if database is None else database
        self.groups = [name.split('_')[-1].split('.')[0] for name in self.fp_in] if groups is None else groups

        # data
        self.data = []
        self.side_chains = []
        self.connections = [(f'[CH3:{self.CONN_MN}]', [0, 3]), (f'[CH3][CH2:{self.CONN_MN}]', [1, 3]),
                            (f'[CH3][CH2][CH2:{self.CONN_MN}]', [2, 3])]

        # flags
        self.json_flag = json_flag
        self.db_flag = db_flag

    @property
    def fp_in(self):
        return self._fp_in

    @property
    def fp_out(self):
        return self._fp_out

    @property
    def mol_db(self):
        return self._mol_db

    @property
    def groups(self):
        return self._groups

    @property
    def data(self):
        return self._data

    @property
    def side_chains(self):
        return self._side_chains

    @property
    def connections(self):
        return self._connections

    @fp_in.setter
    def fp_in(self, val):
        if not isinstance(val, Iterable) and not all(isinstance(fp, str) for fp in val):
            raise TypeError('fp_in must be iterable containing all strings.')
        self._fp_in = val

    @fp_out.setter
    def fp_out(self, val):
        if not isinstance(val, str):
            raise TypeError('fp_out must be a string.')
        self._fp_out = val

    @mol_db.setter
    def mol_db(self, val):
        if not isinstance(val, MolDatabase):
            raise TypeError('mol_db must be an instance of MolDatabase.')
        self._mol_db = val

    @groups.setter
    def groups(self, val):
        if not isinstance(val, Iterable) and not all(isinstance(group, str) for group in val):
            raise TypeError('groups must be an iterable containing all strings.')
        self._groups = val

    @data.setter
    def data(self, val):
        if not isinstance(val, Iterable) and not all(isinstance(d, dict) and d.values() == [
                'heterocycle', 'parent', 'modifications', 'group'] for d in val):
            raise TypeError('data must be an iterable containing dictionaries with proper field names.')
        self._data = val

    @side_chains.setter
    def side_chains(self, val):
        if not isinstance(val, Iterable) and not all(isinstance(mol, Chem.Mol) for mol in val):
            raise TypeError('side_chains must be an iterable containing all rdkit Mols.')
        self._side_chains = val

    @connections.setter
    def connections(self, val):
        if not isinstance(val, Iterable) and not all(isinstance(tup, tuple) and isinstance(tup[0], str) and isinstance(
                tup[1], Iterable) and all(isinstance(key, int) for key in tup[1]) for tup in val):
            raise TypeError('connections must be an iterable containing all tuples with strings and iterable '
                            'of intergers.')
        self._connections = val

    def diversify(self):
        """
        Public instance method. Main driver of class functionality. Calls self.alternate_connection_point() on each
        molecule in self.side_chains with each connection type in self.connections and stores results in self.data.

        Returns:
            bool: True if successful
        """

        for mol, group in self.side_chains:
            for connection, modification in self.connections:
                unique_mols = self.alternate_connection_point(mol, connection)
                self.accumulate_mols(unique_mols, Chem.MolToSmiles(mol), modification, group)

        return True

    def save_data(self):
        """
        Public instance method. Writes data defined in self.data to json file defined by self.fp_out and to 'molecules'
        database in the 'side_chains' collection. Must have corresponding flags set to True (self.json_flag
        and self.db_flag).

        Returns:
            bool: True if successful
        """

        if self.json_flag:
            self.write_json(self.data, self.fp_out)

        if self.db_flag:
            self.mol_db.insert_side_chains(self.data)

        return True

    def get_side_chains(self):
        """
        Reads all side_chains from files in self.fp_in into self.side_chains.

        Returns:
            bool: True if successful
        """

        for fp, group in zip(self.fp_in, self.groups):
            for mol in read_mols(fp):
                self.side_chains.append((mol, group))

        return True

    def alternate_connection_point(self, mol, connection):
        """
        Creates a set of new molecules by attaching an alkyl chain (which becomes the attachment point to peptide
        backbone) to every eligble position on the side chain. Eligiblity of an atom is defined as:
            Carbon - Must have 1 or 2 hydrogens
            Nitrogen, Oxygen, Sulfur - Must have > 0 hydrogens

        Args:
            mol (rdkit Mol): The side chain molecule
            connection (str): The atom mapped smarts SMARTS string of the alkyl attachment chain

        Raises:
            SystemExit: If connection is not atom mapped

        Returns:
            set: A set of unique SMILES strings representing the side chain with alkyl chains attached at different
                positions
        """

        mols = set()
        attached = set()
        connection = Chem.MolFromSmarts(connection)

        # check if connecting atom is atom mapped
        map_nums = [atom.GetAtomMapNum() for atom in connection.GetAtoms()]
        if self.CONN_MN not in map_nums:
            raise SystemExit('Need to specifiy connecting atom with atom map number')

        # try to make attachment at each atom
        for atom in mol.GetAtoms():

            # detetmine atom eligibility
            atom_idx = None
            found = False
            if atom.GetSymbol() == 'C' and atom.GetTotalNumHs() != 0 and atom.GetTotalNumHs() < 3:
                atom.SetAtomMapNum(self.HETERO_MN)
                atom_idx = atom.GetIdx()
                attached.add(atom_idx)
                found = True
            elif atom.GetSymbol() in ('N', 'O', 'S') and atom.GetTotalNumHs() != 0:
                atom.SetAtomMapNum(self.HETERO_MN)
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
                if combo_atom.GetAtomMapNum() == self.HETERO_MN:
                    mol_atom = combo_atom.GetIdx()
                    combo_atom.SetAtomMapNum(0)
                    Chem.Mol.GetAtomWithIdx(mol, atom_idx).SetAtomMapNum(0)
                elif combo_atom.GetAtomMapNum() == self.CONN_MN:
                    conn_atom = combo_atom.GetIdx()
                    combo_atom.SetAtomMapNum(0)

            # create bond
            combo.AddBond(mol_atom, conn_atom, order=Chem.rdchem.BondType.SINGLE)

            # fix hydrogen counts
            atom_react = combo.GetAtomWithIdx(mol_atom)
            if atom_react.GetSymbol() in ('N', 'O', 'S'):
                atom_react.SetNumExplicitHs(0)
            elif atom_react.GetSymbol() == 'C' and Chem.Atom.GetTotalNumHs(atom_react) > 0:
                atom_react.SetNumExplicitHs(Chem.Atom.GetTotalNumHs(atom_react) - 1)

            try:
                Chem.SanitizeMol(combo)
                mols.add(Chem.MolToSmiles(combo))
            except ValueError:
                print('Sanitize Error!')
                print('Side Chain:', Chem.MolToSmiles(mol) + '\n')

        return mols

    def accumulate_mols(self, mols, parent, modifications, group):
        """
        Stores all data associated with the modified side_chain in a dictionary and appends it to self.data.

        Args:
            mols (iterable): A set containing the unique SMILES strings of the modified side_chain
            parent (rdkit Mol): The parent heterocycle from which the modified heterocycle was derived
            modifications (list): A list containing indicators of the modifications that were made to the heterocycle
            group (str): The group (sdf file) the parent heterocycle came from
        """

        for smiles in mols:
            doc = {}
            doc['heterocycle'] = smiles
            doc['parent'] = parent
            doc['modifications'] = modifications
            doc['group'] = group
            self.data.append(doc)


def main():
    """
    Driver function. Parses arguments, constructs class, and performs operations on data.
    """

    parser = argparse.ArgumentParser(description='Creates a unique set of molecules by attaching varying length alkyl '
                                     'chains to all elgible positions on the heterocycle. Alkyl chains include '
                                     "methyl, ethyl, and propyl. The last word in the input files' name dictates which "
                                     'group the heterocycle belongs to.')
    parser.add_argument('--f_in', dest='input', nargs='+',
                        default=['side_chains_likely1.sdf'], help='The sdf file(s) containing monomer side chains.')
    parser.add_argument('--f_out', dest='output', default='side_chains.json', help='The output json file.')
    parser.add_argument('--fp_in', dest='fp_in', default='chemdraw/pre_monomer/',
                        help='The input filepath relative to script')
    parser.add_argument('--fp_out', dest='fp_out', default='smiles/pre_monomer/',
                        help='The ouput filepath relative to script')
    parser.add_argument('--db', dest='database', default='molecules',
                        help='The mongoDB database to connect to')
    parser.add_argument('--host', dest='host', default='localhost',
                        help='The host MongoDB server to connect to')
    parser.add_argument('--port', dest='port', type=int, default=27017,
                        help='The port on host server to connect to')

    args = parser.parse_args()

    # get to input and output files
    fp_in = [join(args.fp_in, file) for file in args.input]
    fp_out = join(args.fp_out, args.output)

    # create classes and perform operations
    mol_db = MolDatabase(database=args.database, host=args.host, port=args.port)
    modifier = SideChainModifier(fp_in, fp_out, database=mol_db)
    if modifier.get_side_chains() and modifier.diversify():
        return modifier.save_data()

    return False


if __name__ == '__main__':
    main()
