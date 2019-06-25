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

from macrocycles.config import (RG_DOC_TYPE, RG_INPUT_COL, RG_INPUT_DIR,
                                RG_OUTPUT_COL, RG_OUTPUT_DIR, RXN_MAP_NUM, TEMPLATE_RXN_MAP_NUM)
from macrocycles.utils import (Base, Flags, IOPaths, MongoParams,
                               create_logger, set_flags)

LOGGER = create_logger(__name__, INFO)


class ReactionGenerator(Base):
    """
    Class for generating atom mapped reaction SMARTS strings from atom mapped template and side chain SMILES strings.
    Inherits from Base.

    Attributes:
        side_chains (list): Contains the different atom mapped side chains.
        templates (list): Contains the different atom mapped templates.
        reaction (str): The type of reaction that is being performed (i.e. Friedel-Crafts, Pictet-Spengler,..).
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
        f_in = [os.path.join(RG_INPUT_DIR, file) for file in kwargs['f_in']] if 'f_in' in kwargs else ['']
        f_out = os.path.join(RG_OUTPUT_DIR, kwargs['f_out']) if 'f_out' in kwargs else ''
        mongo_params = kwargs['mongo_params'] if 'mongo_params' in kwargs else MongoParams(
            RG_INPUT_COL, RG_DOC_TYPE, RG_OUTPUT_COL)
        mongo_params = mongo_params if not no_db else None
        super().__init__(IOPaths(f_in, f_out), mongo_params, LOGGER, input_flags, output_flags)

        # data
        self.side_chains = []
        self.templates = []
        self.reaction = None

    def load_data(self):
        """
        Overloaded method for loading input data.

        Returns:
            bool: True if successful.
        """

        try:
            self.side_chains, self.templates = super().load_data()
            self.templates = list(self.templates)  # convert smaller of two generators to list
        except IndexError:
            self.logger.exception('Check MongoParams contains correct number of input_cols and input_types. '
                                  f'self.mongo_params = {self.mongo_params}')
        except Exception:
            self.logger.exception('Unexpected exception occured.')
        else:
            return True

        return False

    def generate_reactions(self):
        """
        Top level function that generates reaction SMARTS for all pairs of templates and side_chains using helper
        function generate_rxn_template().

        Args:
            templates (iterable): An iterable object containing template documents such as pymongo cursor on the templates
                collection
            side_chains (iterable): An iterable object containing side chain documents such as a pymongo cursor on the
                side_chains colelction
            show_progress (bool, optional): Show progress bar. Defaults to False.

        Returns:
            generator: Generator object of tuples containing the reactions SMARTS string, template document, and side_chain
                document
        """

        try:
            for side_chain in self.side_chains:
                sc_mol = Chem.MolFromSmiles(side_chain['smarts'])
                for template in self.templates:
                    temp_mol = Chem.MolFromSmiles(template['smarts'])
                    reaction = self.friedel_crafts(temp_mol, sc_mol)
                    self.accumulate_data(template, side_chain, reaction)
        except Exception:
            self.logger.exception('Unexpected exception occured.')
        else:
            self.logger.info(f'Successfully generated {len(self.result_data)} reactions')
            return True

        return False

    def friedel_crafts(self, template, side_chain):
        """
        Generate a reaction SMARTS string for Friedel-Crafts reaction with the template and side_chain.

        Args:
            template (rdkit Mol): The template's atom mapped SMILES string
            side_chain (rdkit Mol): The side chain's atom mapped SMILES string

        Returns:
            str: The reaction SMARTS string
        """

        self.reaction = 'friedel_crafts'
        template = Chem.DeleteSubstructs(template, Chem.MolFromSmiles(
            'CC(C)(C)OC(=O)O'))  # t-butyl carbonate leaving group

        try:
            product = Base.merge(template, side_chain, TEMPLATE_RXN_MAP_NUM, RXN_MAP_NUM, clear_map_nums=False)
        except ValueError:
            self.logger.exception(f'Sanitize Error! template = {template}, side_chain = {side_chain}')
        else:
            product = Chem.MolToSmiles(product)
            template = Chem.MolToSmiles(template)
            side_chain = Chem.MolToSmiles(side_chain)
            return '(' + template + '.' + side_chain + ')>>' + product

    def accumulate_data(self, template, side_chain, reaction):
        """
        Stores all data associated with the different reactions into a dictionary and appends it so self.result_data.

        Args:
            template (dict): The associated data of the templates.
            side_chain (dict): The associated data of the side chain that the regioisomers are derived from.
            reaction (str): The atom mapped reaction SMARTS string.
        """

        doc = OrderedDict([('ID', template['ID'] + side_chain['ID']),
                           ('type', 'reaction'),
                           ('smarts', reaction),
                           ('template', template['ID']),
                           ('side_chain', side_chain),
                           ('reaction', self.reaction)])
        self.result_data.append(doc)


def main():
    """
    Driver function. Parses arguments, constructs class, and performs operations on data.
    """

    parser = argparse.ArgumentParser(description='Generates and stores a SMARTS reaction template based on template '
                                     'and side chain SMILES strings stored in the MongoDB database. The product in the '
                                     'reaction template corresponds to the connection of the template atom with atom '
                                     'map number = 1 to the side chain atom with atom map number = 2. The SMARTS '
                                     'reaction template can then be applied to template-peptide linked molecules to '
                                     'form a macrocycle.')
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
    generator = ReactionGenerator(f_in=f_in, f_out=args.f_out, input_flags=input_flags,
                                  output_flags=output_flags, no_db=args.no_db)
    if generator.load_data() and generator.generate_reactions():
        # return generator.save_data()
        print(len(generator.result_data))
        print(generator.result_data[0])

    return False


if __name__ == '__main__':
    main()
