import argparse
import json
from pathlib import Path
from copy import deepcopy
from logging import INFO
import os
import sys
from pprint import pprint

from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from tqdm import tqdm

from macrocycles.utils import Base, IOPaths, Flags, create_logger, MongoParams, set_flags
from macrocycles.config import CE_INPUT_DIR, CE_OUTPUT_DIR, CE_INPUT_COL, CE_DOC_TYPE, CE_OUTPUT_COL

# TODO: Can make enumeration process quicker by only searching for reaction templates that correspond to side chains in
# the peptide. Use monomer information in template_peptide documents to traceback side_chain structure and query db to
# only find reaction templates that contain that side chain.

LOGGER = create_logger(__name__, INFO)


class CandidateEnumerator(Base):
    """
    Class for enumerating the candidate macrocycles by applying reactions to the tp_hybrids. Inherits from Base.

    Attributes:
        tp_hybrids (list): Contains the different side chains for generating regioisomers.
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
        f_in = [os.path.join(CE_INPUT_DIR, file) for file in kwargs['f_in']] if 'f_in' in kwargs else ['']
        f_out = os.path.join(CE_OUTPUT_DIR, kwargs['f_out']) if 'f_out' in kwargs else ''
        mongo_params = kwargs['mongo_params'] if 'mongo_params' in kwargs else MongoParams(
            CE_INPUT_COL, CE_DOC_TYPE, CE_OUTPUT_COL)
        mongo_params = mongo_params if not no_db else None
        super().__init__(IOPaths(f_in, f_out), mongo_params, LOGGER, input_flags, output_flags)

        # data
        self.tp_hybrids = []

    def load_data(self):
        """
        Overloaded method for loading input data.

        Returns:
            bool: True if successful.
        """

        try:
            self.tp_hybrids = super().load_data()[0]
        except IndexError:
            self.logger.exception('Check MongoParams contains correct number of input_cols and input_types. '
                                  f'self.mongo_params = {self.mongo_params}')
        except Exception:
            self.logger.exception('Unexpected exception occured.')
        else:
            return True

        return False

    def generate_candidates(self):

        for tp_hybrid in self.tp_hybrids:
            reactant = tp_hybrid['smiles']
            print('', reactant)
            for i, monomer in enumerate(set(tp_hybrid['peptide']['monomers'])):
                print('\t', monomer)
                monomer = self.mongo_db['molecules'].find_one({'ID': monomer})
                # print(monomer)
                reactions = list(self.mongo_db['reactions'].find(
                    {'type': 'reaction', 'side_chain.side_chain.parent.ID': monomer['side_chain']['parent']['ID']}))
                reactions.extend(list(self.mongo_db['reactions'].find(
                    {'type': 'reaction', 'side_chain': monomer['side_chain']['ID']}
                )))
                # print(reactions[0])
                print(len(reactions))
                # print(reactions)
                products, reacting_side_chains, atom_idxs = self.apply_reaction(reactant, reactions)
                self.accumulate_data(products, reacting_side_chains, atom_idxs, tp_hybrid, reactions, i)

        return True

    def apply_reaction(self, reactant, reactions):
        """
        Applys all reaction templates to a reactant and collects all unique products

        Args:
            reactant (str): The reactant SMILES string
            rxn_docs (pymongo cursor): The collection of reaction documents containing reaction template SMARTS strings

        Returns:
            list: A list of SMARTS strings of all unique products
            list: A list of side chains that reacted in order to form the products
            list: A list of reacting atom indices in the side_chains that formed products
        """

        # print(reactant)
        reactant = Chem.MolFromSmiles(reactant)

        # apply reactions
        reaction_info = {}
        for doc in reactions:
            rxn = AllChem.ReactionFromSmarts(doc['smarts'])
            prod = rxn.RunReactants((reactant,))
            # print(prod)

            # get unique products
            for mol in prod:
                Chem.SanitizeMol(mol[0])
                smiles = Chem.MolToSmiles(mol[0])
                # print(smiles)
                if smiles not in reaction_info.keys():
                    try:
                        reaction_info[smiles] = (doc['side_chain'], doc['side_chain']['rxn_atom_idx'])
                    except:
                        reaction_info[smiles] = (doc['side_chain'], None)

        # print(reaction_info)
        unique_prods, sc_and_idx = zip(*[(key, value) for key, value in reaction_info.items()])
        reacting_side_chains, atom_idxs = zip(*list(sc_and_idx))
        # print(reacting_side_chains)
        # print(atom_idxs)

        return unique_prods, reacting_side_chains, atom_idxs

    def accumulate_data(self, macrocycles, reacting_side_chains, atom_idxs, tp_hybrid, reactions, i):

        # for j, (macrocycle, reacting_side_chain, atom_idx) in enumerate(zip(macrocycles, reacting_side_chains, atom_idxs)):
        doc = {
            'ID': 'c' + str(i) + tp_hybrid['ID'],
            'type': 'candidate',
            'smiles': macrocycles,
            'kekule': '',
            'reacting_side_chains': reacting_side_chains,
            'atom_idxs': atom_idxs,
            'tp_hybrid': tp_hybrid
            # 'reaction': reactions
        }
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
    enumerator = CandidateEnumerator(f_in=f_in, f_out=args.f_out, input_flags=input_flags,
                                     output_flags=output_flags, no_db=args.no_db)
    if enumerator.load_data() and enumerator.generate_candidates():
        return enumerator.save_data()
        # print(len(enumerator.result_data))
        # pprint(enumerator.result_data[-2:])

    return False

    # parser = argparse.ArgumentParser(description='')
    # parser.add_argument('-f', '--fin', dest='in_file', nargs='+', default=['length3_temp1.json'],
    #                     help='The input json file(s) containing the SMILES strings of template-peptide connected molecule')
    # parser.add_argument('-fp', '--fp_in', dest='fp_in', default='smiles/template_peptide/c_term/',
    #                     help='The filepath to the input files relative to base project directory')
    # parser.add_argument('-d', '--db', dest='database', default='rxn_templates',
    #                     help='The mongoDB database to connect to')
    # parser.add_argument('-hn', '--host', dest='host', default='localhost',
    #                     help='The host MongoDB server to connect to')
    # parser.add_argument('-p', '--port', dest='port', type=int, default=27017,
    #                     help='The port on host server to connect to')
    # parser.add_argument('-s', '--store', dest='store', action='store_false', help='Toggle to not store results')
    # parser.add_argument('--show_progress', dest='progress', action='store_false',
    #                     help='Show progress bar. Defaults to False')

    # args = parser.parse_args()

    # # generate filepath to the input file
    # base_path = Path(__file__).resolve().parents[1]
    # fp_in = [str(base_path / args.fp_in / file) for file in args.in_file]

    # # establish database connection and retrieve reaction documents
    # db = Database(host=args.host, port=args.port, db=args.database)
    # rxn_docs = db.find_all('reactions')

    # # for each reactant generate candidates
    # db.db = db.client['molecules']
    # for fp in fp_in:
    #     with open(fp, 'r') as f:
    #         for doc in tqdm(json.load(f), disable=args.progress):

    #             # extract data
    #             reactant = doc['template-peptide']
    #             peptide = doc['peptide']
    #             template = doc['template']
    #             monomers = doc['monomers']

    #             products, reacting_side_chains, atom_idx = generate_candidates(reactant, deepcopy(rxn_docs))

    #             if args.store:
    #                 db.insert_candidates(reactant, products, len(products), peptide,
    #                                      template, monomers, reacting_side_chains, atom_idx)


if __name__ == '__main__':
    main()
