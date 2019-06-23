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

from macrocycles.config import (BB_MAP_NUM, MG_DOC_TYPE, MG_INPUT_COL,
                                MG_INPUT_DIR, MG_OUTPUT_COL, MG_OUTPUT_DIR,
                                SC_MAP_NUM, STEREOCHEMISTRY)
from macrocycles.exceptions import AtomSearchError
from macrocycles.utils import (Base, Flags, IOPaths, MongoParams, Smiles,
                               create_logger, set_flags)

LOGGER = create_logger(name=__name__, level=INFO)


class MonomerGenerator(Base):
    """
    Class for combining side_chains and amino acid backbones to form monomers. Inherits from Base.

    Attributes:
        backbones (list): Contains the different amino acid backbones.
        side_chains (list): Contains the side_chains and the associated data as dictionaries.
    """

    def __init__(self, required, logger=LOGGER, input_flags=Flags(False, False, True, False),
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
        f_in = [os.path.join(MG_INPUT_DIR, file) for file in kwargs['f_in']] if 'f_in' in kwargs else ['']
        f_out = os.path.join(MG_OUTPUT_DIR, kwargs['f_out']) if 'f_out' in kwargs else ''
        mongo_params = kwargs['mongo_params'] if 'mongo_params' in kwargs else MongoParams(
            MG_INPUT_COL, MG_DOC_TYPE, MG_OUTPUT_COL)
        mongo_params = mongo_params if not no_db else None
        super().__init__(IOPaths(f_in, f_out), mongo_params, LOGGER, input_flags, output_flags)

        # data
        self.backbones = []
        self.side_chains = []
        self.stereo = STEREOCHEMISTRY
        self.required = required

    def load_data(self):

        try:
            self.side_chains, self.backbones = super().load_data()
        except ValueError:
            self.logger.exception('Check MongoParams contains correct number of input_cols and input_types. '
                                  f'self.mongo_params = {self.mongo_params}')
            return False
        except Exception:
            self.logger.exception('Unexcepted exception occured.')
            return False
        else:
            return True

    def generate_monomers(self):
        """
        Top level function that calls self.create_monomer() on each side_chain - backbone - stereochemistry combination
        and calls self.accumulate_mols() with the resulting monomers.

        Returns:
            bool: True if successful.
        """

        try:
            for backbone in self.backbones:
                mol_bb = Chem.MolFromSmarts(backbone['smiles'])
                for side_chain in self.side_chains:
                    monomers = []
                    mol_sc = Chem.MolFromSmiles(side_chain['smiles'])
                    for stereo in self.stereo:
                        monomers.extend(self.create_monomer(mol_bb, mol_sc, stereo))
                    self.accumulate_mols(monomers, backbone, side_chain)
        except AtomSearchError:
            self.logger.exception()
            return False
        else:
            return True

    def create_monomer(self, backbone, side_chain, stereo):
        """
        Connects the side chain to the backbone at the designated attachment points.

        Args:
            backbone (rdkit Mol): The backbone to which the side chain will be attached to.
            side_chain (rdkit Mol): The side chain being attached to the backbone structure.
            stereo (str): The stereochemistry at the attachment point of the backbone and side chain.

        Returns:
            rdkit Mol: The resulting monomer from connecting the backbone and side chain molecules
        """

        # attachment point at terminal end of alkyl chain
        patt = Chem.MolFromSmarts('[CH3]')
        matches = sorted(side_chain.GetSubstructMatches(patt, useChirality=False), key=lambda x: x[0], reverse=True)

        if not len(matches):
            raise AtomSearchError(f'Couldn\'t find attchment point of: {Chem.MolToSmiles(side_chain)}')

        monomers = []
        for pairs in matches:

            # set atom map number for attachment point of side chain
            sc_atom_idx = None
            for atom_idx in pairs:
                atom = Chem.Mol.GetAtomWithIdx(side_chain, atom_idx)
                if atom.GetSymbol() == 'C' and Chem.Atom.GetTotalNumHs(atom) == 3 and not Chem.Atom.IsInRing(atom):
                    atom.SetAtomMapNum(SC_MAP_NUM)
                    sc_atom_idx = atom_idx
                    break
            else:
                raise AtomSearchError(f'Couldn\'t find attchment point of: {Chem.MolToSmiles(side_chain)}')

            combo = Chem.RWMol(Chem.CombineMols(side_chain, backbone))

            # get reacting atom indicies
            bb_atom = None
            sc_atom = None
            for atom in combo.GetAtoms():
                if atom.GetAtomMapNum() == BB_MAP_NUM:
                    bb_atom = atom.GetIdx()
                    atom.SetAtomMapNum(0)
                elif atom.GetAtomMapNum() == SC_MAP_NUM:
                    sc_atom = atom.GetIdx()
                    atom.SetAtomMapNum(0)
                    Chem.Mol.GetAtomWithIdx(side_chain, sc_atom_idx).SetAtomMapNum(0)

            combo.AddBond(bb_atom, sc_atom, order=Chem.rdchem.BondType.SINGLE)

            # adjust Stereochemistry, fix hydrogen counts
            atom_react = combo.GetAtomWithIdx(bb_atom)
            atom_react.SetNumExplicitHs(Chem.Atom.GetTotalNumHs(atom_react) - 1)
            if stereo == 'CCW':
                atom_react.SetChiralTag(Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW)
            else:
                atom_react.SetChiralTag(Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW)

            try:
                Chem.SanitizeMol(combo)
                Chem.Kekulize(combo)
            except ValueError:
                self.logger.exception(
                    f'Sanitize Error! Side Chain: {Chem.MolToSmiles(side_chain)}, backbone: {backbone}')
            else:
                monomers.append(Smiles(Chem.MolToSmiles(combo), Chem.MolToSmiles(combo, kekuleSmiles=True)))
                self.logger.debug(f'Success! Monomer {Chem.MolToSmiles(combo)}')
                attach_complete = True

        return monomers

    def accumulate_mols(self, monomers, backbone, side_chain):
        """
        Stores all data associated with the monomer in a dictionary and appends it to self.result_data.

        Args:
            monomers (iterable): Contains Smiles namedtuples that contain the SMILES and kekule SMILES of the monomers.
            backbone (dict): The dictionary containing the associated backbone data.
            side_chain (dict): The dictionary containing the associated side chain data.
        """

        for i, monomer in enumerate(monomers):
            doc = OrderedDict([('ID', 'm' + str(i) + side_chain['ID']),
                               ('type', 'monomer'),
                               ('smiles', monomer.smiles),
                               ('kekule', monomer.kekule),
                               ('backbone', backbone['ID']),
                               ('side_chain', side_chain),
                               ('required', self.required),
                               ('group', side_chain['parent']['group'])
                               ])
            self.result_data.append(doc)


def main():
    """
    Driver function. Parses arguments, constructs class, and performs operations on data.
    """

    parser = argparse.ArgumentParser(description='Generates monomers by combining side chains and backbone together. If '
                                     'the set of generated monomers are to be required in peptide generation, need to '
                                     'spcify option --required.')
    parser.add_argument('input', choices=['json', 'txt', 'mongo', 'sql'],
                        help='Specifies the format that the input data is in.')
    parser.add_argument('output', nargs='+', choices=['json', 'txt', 'mongo', 'sql'],
                        help='Specifies what format to output the result data in.')
    parser.add_argument('--bb_file', dest='bb_file', required='json' in sys.argv or 'txt' in sys.argv,
                        help='The input file containing backbone data relative to default input directory defined '
                        'in config.py.')
    parser.add_argument('--sc_file', dest='sc_file', required='json' in sys.argv or 'txt' in sys.argv,
                        help='The input file containing side chain data relative to default input directory defined '
                        'in config.py.')
    parser.add_argument('--f_out', dest='f_out', default='custom',
                        help='The output file relative to the default output directory defined in config.py.')
    parser.add_argument('--no_db', dest='no_db', action='store_true',
                        help='Turns off default connection that is made to the database.')
    parser.add_argument('--show_progress', dest='progress', action='store_false',
                        help='Show progress bar. Defaults to False')
    parser.add_argument('--required', dest='required', action='store_true',
                        help='Determines whether the monomer will be in the required pool for generating peptides. '
                        'Defaults to False.')

    args = parser.parse_args()

    # check for proper file specifications
    if args.input in ['json', 'txt']:
        extensions = [args.bb_file.split('.')[-1]]
        extensions.append(args.sc_file.split('.')[-1])
        if not all([args.input == extension for extension in extensions]):
            LOGGER.error('File extension of the input files does not match the specified format')
            raise OSError('File extension of the input files does not match the specified format')

    # configure I/O
    input_flags, output_flags = set_flags(args.input, args.output)
    f_in = [args.sc_file, args.bb_file] if args.input in ['json', 'txt'] else ['']

    generator = MonomerGenerator(required=args.required, f_in=f_in, f_out=args.f_out, input_flags=input_flags,
                                 output_flags=output_flags, no_db=args.no_db)
    if generator.load_data() and generator.generate_monomers():
        return generator.save_data()

    return False


if __name__ == '__main__':
    main()
