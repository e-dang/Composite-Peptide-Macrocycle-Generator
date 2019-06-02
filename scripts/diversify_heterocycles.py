"""
Written by Eric Dang.
github: https://github.com/e-dang
email: erickdang@g.ucla.edu
"""

import argparse

from rdkit import Chem

from base import Base
from database import MolDatabase
from utils import read_mols

# TODO: allow for side chains with methyl group already attached, such as 5-methyl indole.


class Diversify_Heterocycles(Base):
    """
    Class for attaching varying length alkyl chains / attachment points to heterocycles.

    Attributes:
        fp_in (list): Contains all filepaths to input files containing the parent heterocycles.
        fp_out (str): The filepath to the output json file.
        mol_db (MolDatabase): An instance of MolDatabase.
        groups (list): Contains the groups from which each set of heterocycles comes from.
        connections (list): Contains tuples, each of which contain atom mapped SMARTS strings of the alkyl attachment
            chain and the corresponding modification array.
        data (list): Contains the dictionaries associated with the modified heterocycles.
        json_flag (bool): Determines whether to write data to the json file.
        db_flag (bool): Determines whether to write data to the database.
    """

    HETERO_MN = 1
    CONN_MN = 2

    def __init__(self, fp_in, fp_out, groups, database='molecules', host='localhost', port=27017, json_flag=True,
                 db_flag=True):
        """
        Constructor.

        Args:
            fp_in (list): Contains all filepaths to input files containing the parent heterocycles.
            fp_out (str): The filepath to the output json file.
            groups (list): Contains the groups from which each set of heterocycles comes from.
            database (str, optional): The database to connect to. Defaults to 'molecules'.
            host (str, optional): The host containing the database to connect to. Defaults to 'localhost'.
            port (int, optional): The port to connect to. Defaults to 27017.
            json_flag (bool): Determines whether to write data to the json file.
            db_flag (bool): Determines whether to write data to the database.
        """

        super().__init__()
        self.fp_in = [str(self.project_dir / file) for file in fp_in]
        self.fp_out = str(self.project_dir / fp_out)
        self.mol_db = MolDatabase(database=database, host=host, port=port)
        self.groups = groups
        self.connections = [(f'[CH3:{self.CONN_MN}]', [0, 3]), (f'[CH3][CH2:{self.CONN_MN}]', [1, 3]),
                            (f'[CH3][CH2][CH2:{self.CONN_MN}]', [2, 3])]
        self.data = []
        self.json_flag = json_flag
        self.db_flag = db_flag

    def diversify_heterocycles(self):
        """
        Public instance method. Main driver of class functionality. Reads in molecules in self.fp_in and calls
        __alternate_connection_point() on each one with each connection type and stores results in self.data.

        Returns:
            bool: True if successful
        """

        for fp, group in zip(self.fp_in, self.groups):
            for mol in read_mols(fp):
                for connection, modification in self.connections:
                    unique_mols = self.__alternate_connection_point(mol, connection)
                    self.__accumulate_mols(unique_mols, Chem.MolToSmiles(mol), modification, group)

        return True

    def save_data(self):
        """
        Public instance method. Writes data defined in self.data to json file defined by self.fp_out and to 'molecules'
        database in the 'heterocycles' collection. Must have corresponding flags set to True (self.json_flag
        and self.db_flag).

        Returns:
            bool: True if successful
        """

        if self.json_flag:
            self.write_json(self.data, self.fp_out)

        if self.db_flag:
            self.mol_db.insert_heterocycles(self.data)

        return True

    def __alternate_connection_point(self, mol, connection):
        """
        Private instance method. Creates a set of new molecules by attaching an alkyl chain (which becomes the attachment
        point to peptide backbone) to every eligble position on the side chain. Eligiblity of an atom is defined as:
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
            print('Need to specifiy connecting atom with atom map number')
            raise SystemExit

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

    def __accumulate_mols(self, mols, parent, modifications, group):
        """
        Private instance method. Stores all data associated with the modified heterocycles in a dictionary and appends
        it to self.data.

        Args:
            mols (iterable): A set containing the unique SMILES strings of the modified heterocycles
            parent (rdkit Mol): The parent heterocycle from which the modified heterocycle was derived
            modifications (list): A list containing indicators of the modifications that were made to the heterocycle
            group (str): The group (sdf file) the parent heterocycle came from
        """

        for smiles in mols:
            doc = {}
            doc['heterocycle'] = smiles
            doc['parent'] = parent
            doc['modification'] = modifications
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
    parser.add_argument('--f_out', dest='out', default='side_chains.json', help='The output json file.')
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

    # get side_chain group name from input file name
    groups = [name.split('_')[-1].split('.')[0] for name in args.input]

    # get absolute filepath to input and output files
    fp_in = [args.fp_in + '/' + file for file in args.input]
    fp_out = args.fp_out + '/' + args.out

    # create class and perform operations
    diversifier = Diversify_Heterocycles(fp_in, fp_out, groups, args.database, args.host, args.port)
    if diversifier.diversify_heterocycles():
        return diversifier.save_data()

    return False


if __name__ == '__main__':
    main()
