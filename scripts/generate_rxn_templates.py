import argparse

from rdkit import Chem
from tqdm import tqdm

from utils import Database


def generate_rxn_templates(templates, side_chains, show_progress=False):
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

    # create SMARTS reaction template between each template and side chain
    for template in tqdm(templates, desc='Templates: ', disable=show_progress):
        for side_chain in tqdm(side_chains, desc='Side chains: ', disable=show_progress):
            yield generate_rxn_template(template, side_chain), template, side_chain


def generate_rxn_template(template, side_chain):
    """
    Uses helper function merge() to generate product SMILES string, then use helper function generate_rxn_smarts() to
    create reaction SMARTS string

    Args:
        template (pymongo doc): The template's database document
        side_chain (pumongo doc): The side chain's database document

    Returns:
        str: The full reaction SMARTS string
    """

    # combine SMILES strings to get full reaction template SMARTS string
    return generate_rxn_smarts(template['atom_mapped_smiles'], side_chain['atom_mapped_smiles'],
                               merge(template, side_chain))


def generate_rxn_smarts(template, side_chain, product):
    """
    Generate a reaction SMARTS string from the arugments

    Args:
        template (str): The template's atom mapped SMILES string
        side_chain (str): The side chain's atom mapped SMILES string
        product (str): The product's atom mapped SMILES string

    Returns:
        str: The reaction SMARTS string
    """

    return '(' + template + '.' + side_chain + ')>>' + product


def merge(template, side_chain):
    """
    Connects the template and side chain together at the designated reacting site (marked by atom map numbers), and
    converts the product to a SMILES string.

    Args:
        template (dict): The template's database document
        side_chain (dict): The side chain's database document

    Returns:
        str: SMILES string representation of the product of the template reacting with the side chain
    """

    # convert SMILES strings to mols
    temp = Chem.MolFromSmiles(template['atom_mapped_smiles'])
    sc = Chem.MolFromSmiles(side_chain['atom_mapped_smiles'])

    # remove substruct from template and combine mols
    temp = Chem.DeleteSubstructs(temp, Chem.MolFromSmiles(template['substruct']))
    combo = Chem.RWMol(Chem.CombineMols(temp, sc))

    # get reacting atom indicies
    temp_atom = None
    sc_atom = None
    for atom in combo.GetAtoms():
        if atom.GetAtomMapNum() == int(template['rxn_map_num']):
            temp_atom = atom.GetIdx()
        elif atom.GetAtomMapNum() == int(side_chain['rxn_map_num']):
            sc_atom = atom.GetIdx()

    # fix hydrogen counts
    atom_react = combo.GetAtomWithIdx(sc_atom)
    if atom_react.GetSymbol() == 'N' or atom_react.GetSymbol() == 'O' or atom_react.GetSymbol() == 'S':
        atom_react.SetNumExplicitHs(0)
    elif atom_react.GetSymbol() == 'C' and Chem.Atom.GetTotalNumHs(atom_react) > 0:
        atom_react.SetNumExplicitHs(Chem.Atom.GetTotalNumHs(atom_react) - 1)

    # create bond
    combo.AddBond(temp_atom, sc_atom, order=Chem.rdchem.BondType.SINGLE)

    try:
        Chem.SanitizeMol(combo)
    except ValueError:
        print('Error!')
        print('Template:', template['smiles'])
        print('Side Chain:', side_chain['smiles'])

    return Chem.MolToSmiles(combo)


def main():
    parser = argparse.ArgumentParser(description='Generates and stores a SMARTS reaction template based on template '
                                     'and side chain SMILES strings stored in the MongoDB database. The product in the '
                                     'reaction template corresponds to the connection of the template atom with atom '
                                     'map number = 1 to the side chain atom with atom map number = 2. The SMARTS '
                                     'reaction template can then be applied to template-peptide linked molecules to '
                                     'form a macrocycle.')
    parser.add_argument('templates', nargs='+', choices=[1, 2, 3], type=int,
                        help='Which template(s) SMILES strings to import; 1 - temp1a, 2 - temp2, 3 - temp3')
    parser.add_argument('--db', dest='database', default='rxn_templates',
                        help='The mongoDB database to connect to')
    parser.add_argument('--host', dest='host', default='localhost',
                        help='The host MongoDB server to connect to')
    parser.add_argument('--port', dest='port', type=int, default=27017,
                        help='The port on host server to connect to')
    parser.add_argument('--store', dest='store', action='store_true', help='Store results (defaults to False)')
    parser.add_argument('--silence', dest='silence', action='store_false',
                        help='Disable printing output to console (defaults to False)')
    parser.add_argument('--show_progress', dest='progress', action='store_false',
                        help='Show progress bar (defaults to False)')

    args = parser.parse_args()

    # get template names for db querying
    templates = [name for ind, name in enumerate(['temp1', 'temp2', 'temp3']) if ind + 1 in args.templates]

    # establish connection and retrieve template and side chain data
    db = Database(host=args.host, port=args.port, db=args.database)
    temp_docs = db.find('templates', {'name': {'$in': templates}})
    sc_docs = db.find_all('side_chains')

    for rxn, template, side_chain in generate_rxn_templates(temp_docs, sc_docs, args.progress):

        if args.silence:
            print('Reaction: ', rxn)
            print('Template: ', template['smiles'])
            print('Side Chain: ', side_chain['smiles'])

        if args.store:
            mod_rxn = rxn.replace('\\', '\\\\')
            db.insert_reaction(mod_rxn, template['smiles'], side_chain['smiles'], side_chain['atom_idx'])

    # # create SMARTS reaction template between each template and side chain
    # for temp in tqdm(temp_docs, desc='Templates: ', disable=args.progress):
    #     for sc in tqdm(sc_docs, desc='Side chains: ', disable=args.progress):
    #         rxn = generate_rxn_temp(temp, sc)

    #         if args.store:
    #             new_rxn = rxn.replace('\\', '\\\\')
    #             db.insert_reaction(new_rxn, temp['smiles'], sc['smiles'], sc['atom_idx'])


if __name__ == '__main__':
    main()
