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
from rdkit.Chem import AllChem

from macrocycles.config import (CARBONYL, PEP_MAP_NUM,
                                SUCCINIMIDE, TEMP_MAP_NUM, TP_DOC_TYPE,
                                TP_INPUT_COL, TP_INPUT_DIR, TP_OUTPUT_COL,
                                TP_OUTPUT_DIR)
from macrocycles.utils import (Base, Flags, IOPaths, MongoParams,
                               create_logger, set_flags)


LOGGER = create_logger(__name__, INFO)


class TPHybridGenerator(Base):
    """
    Class for joining templates and peptides together. Inherits from Base.

    Attributes:
        templates (list): Contains the templates and the associated data.
        peptides (list): Contains the peptides and the associated data.
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
        f_in = [os.path.join(TP_INPUT_DIR, file) for file in kwargs['f_in']] if 'f_in' in kwargs else ['']
        f_out = os.path.join(TP_OUTPUT_DIR, kwargs['f_out']) if 'f_out' in kwargs else ''
        mongo_params = kwargs['mongo_params'] if 'mongo_params' in kwargs else MongoParams(
            TP_INPUT_COL, TP_DOC_TYPE, TP_OUTPUT_COL)
        mongo_params = mongo_params if not no_db else None
        super().__init__(IOPaths(f_in, f_out), mongo_params, LOGGER, input_flags, output_flags)

        # data
        self.templates = []
        self.peptides = []

    def load_data(self):
        """
        Overloaded method for loading input data.

        Returns:
            bool: True if successful.
        """

        try:
            self.peptides, self.templates = super().load_data()
            self.templates = list(self.templates)  # convert smaller of two generators to list
        except ValueError:
            self.logger.exception('Check MongoParams contains correct number of input_cols and input_types. '
                                  f'self.mongo_params = {self.mongo_params}')
        except Exception:
            self.logger.exception('Unexcepted exception occured.')
        else:
            return True

        return False

    def modify_template(self, template):
        """
        Helper function of get_template(), which performs a deletion of the substructure corresponding to the leaving
        group upon amide bond formation and the assignment of an atom map number to reacting carbon atom

        Args:
            template (rdkit Mol): The rdkit Mol representation of the template

        Returns:
            rdkit Mol: The modified template
        """

        # set atom map number for connection point between peptide and template
        matches = template.GetSubstructMatches(Chem.MolFromSmarts(SUCCINIMIDE))
        for pairs in matches:
            for atom_idx in pairs:
                atom = template.GetAtomWithIdx(atom_idx)
                if atom.GetSymbol() == 'O':
                    neighbors = [neighbor for neighbor in atom.GetNeighbors()]
                    if len(neighbors) == 2:
                        carbon = neighbors[0] if neighbors[0].GetSymbol() == 'C' else neighbors[1]
                        carbon.SetAtomMapNum(TEMP_MAP_NUM)
                        break

        # remove leaving group substruct
        template = AllChem.DeleteSubstructs(template, Chem.MolFromSmiles(SUCCINIMIDE))

        return template

    def generate_tp_hybrids(self):
        """
        Main driver function that calls self.modify_templates() on each template, and self.connect_template_peptide() on
        each unique peptide-template pair. The resulting tp_hybrid is then passed to self.accumulate_data().

        Returns:
            bool: True if successful.
        """

        templates = [self.modify_template(Chem.MolFromSmiles(template['smiles'])) for template in self.templates]

        for peptide in self.peptides:
            peptide_mol = Chem.MolFromSmiles(peptide['smiles'])
            for template, template_mol in zip(self.templates, templates):
                tp_hybrid = self.connect_template_peptide(template_mol, peptide_mol, peptide['N_term'])
                self.accumulate_data(tp_hybrid, template, peptide)

        return True

    def connect_template_peptide(self, template, peptide, n_term):
        """
        Connects the peptide to the template through amide linkage with any primary amine on the peptide. This means a
        with a lysine will form two tp_hybrids, where one is connected through the N-terminal amine and the amine on the
        lysine.

        Args:
            template (rdkit Mol): The template molecule.
            peptide (rdkit Mol): The peptide molecule.
            n_term (str): The type of backbone the N-terminal monomer is made of.

        Returns:
            dict: A dictionary containing the tp_hybrid SMILES strings as keys and the kekule SMILES strings as values.
        """

        # type of backbone on N-terminus
        if n_term == 'alpha':
            patt = Chem.MolFromSmarts('[NH;R]CC(=O)')
        elif n_term in ('beta2', 'beta3'):
            patt = Chem.MolFromSmarts('[NH;R]CCC(=O)')

        # find any primary amine or proline n-terminus
        matches = list(peptide.GetSubstructMatches(Chem.MolFromSmarts(
            '[$([NH2]);!$([NH2]C(=O)*)]'), useChirality=False))
        matches.extend(list(peptide.GetSubstructMatches(patt, useChirality=False)))

        # for each eligible nitrogen form a connection
        results = {}
        for pairs in matches:

            # assign atom map number
            for atom_idx in pairs:
                atom = Chem.Mol.GetAtomWithIdx(peptide, atom_idx)
                if atom.GetSymbol() == 'N':  # primary amines only
                    atom.SetAtomMapNum(PEP_MAP_NUM)
                    break

            # combine and record results
            try:
                tp_hybrid = Base.merge(template, peptide, TEMP_MAP_NUM, PEP_MAP_NUM)
            except ValueError:
                self.logger.exception(f'Sanitize Error! template = {Chem.MolToSmiles(template)}, '
                                      f'peptide = {Chem.MolToSmiles(peptide)}')
            else:
                atom.SetAtomMapNum(0)
                smiles = Chem.MolToSmiles(tp_hybrid)
                Chem.Kekulize(tp_hybrid)
                results[smiles] = Chem.MolToSmiles(tp_hybrid, kekuleSmiles=True)

        return results

    def accumulate_data(self, tp_hybrid, template, peptide):
        """
        Stores all data associated with the tp_hybrid in a dictionary and appends it to self.result_data.

        Args:
            tp_hybrid (dict): The dictionary containing the tp_hybrid SMILES string as keys and kekule SMILES as values.
            template (dict): The dictionary containing the associated template data.
            peptide (dict): The dictionary containing the associated peptide data.
        """

        for i, (smiles, kekule) in enumerate(tp_hybrid.items()):
            doc = OrderedDict([('ID', template['ID'] + str(i) + peptide['ID']),
                               ('type', 'tp_hybrid'),
                               ('smiles', smiles),
                               ('kekule', kekule),
                               ('peptide', peptide),
                               ('template', {'ID': template['ID']})])
            self.result_data.append(doc)


def main():
    """
    Driver function. Parses arguments, constructs class, and performs operations on data.
    """

    parser = argparse.ArgumentParser(description='Connects each peptide defined in the input file to a template and '
                                     'write the resulting molecule to file as a SMILES string')
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

    generator = TPHybridGenerator(f_in=f_in, f_out=args.f_out, input_flags=input_flags,
                                  output_flags=output_flags, no_db=args.no_db)
    if generator.load_data() and generator.generate_tp_hybrids():
        # print(len(generator.result_data))
        # print(generator.result_data[:3])
        return generator.save_data()

    return False


if __name__ == '__main__':
    main()
