from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from utils import Database

def merge(template, side_chain):
    """
    Takes template and side_chain mongo docs, merges the SMARTS strings together at the designated reacting site, and returns the merged SMARTS string.
    """

    # convert smarts strings to mols
    temp = Chem.MolFromSmiles(template['smarts'])
    sc = Chem.MolFromSmiles(side_chain['smarts'])


    # get atom map number of reacting site
    temp_map_num = int(template['rxn_ind'])
    sc_map_num = int(side_chain['rxn_ind'])

    # remove substruct and combine mols
    temp = Chem.DeleteSubstructs(temp, Chem.MolFromSmiles(template['substruct']))
    combo = Chem.RWMol(Chem.CombineMols(temp, sc))

    # get reacting atom indicies
    temp_atom = None
    sc_atom = None
    for atom in combo.GetAtoms():
        if atom.GetAtomMapNum() == temp_map_num:
            temp_atom = atom.GetIdx()
        elif atom.GetAtomMapNum() == sc_map_num:
            sc_atom = atom.GetIdx()

    # check if reacting atom is a nitrogen and if so remove all hydrogens
    atom_react = combo.GetAtomWithIdx(sc_atom)
    if atom_react.GetSymbol() == 'N' or atom_react.GetSymbol() == 'O':
        atom_react.SetNumExplicitHs(0)
    elif atom_react.GetSymbol() == 'C' and Chem.Atom.GetTotalNumHs(atom_react) > 0:
        atom_react.SetNumExplicitHs(Chem.Atom.GetTotalNumHs(atom_react) - 1)


    # create bond
    combo.AddBond(temp_atom, sc_atom, order=Chem.rdchem.BondType.SINGLE)
    Chem.SanitizeMol(combo)

    return Chem.MolToSmiles(combo)

def generate_rxn_temp(template, side_chain, store=False, verbose=False):
    """
    Takes template and side_chain mongo docs and combines them to create and return a SMARTS reaction template.
    """

    # get smarts strings
    temp = template['smarts']
    sc = side_chain['smarts']
    prod = merge(template, side_chain)

    # combine smarts strings
    rxn = '(' + temp + '.' + sc + ')>>' + prod

    if verbose:
        Draw.ReactionToImage(AllChem.ReactionFromSmarts(rxn), subImgSize=(500,500)).show()

    if store:
        new_rxn = rxn.replace('\\', '\\\\')
        db = Database()
        db.insert_reaction(template['name'], side_chain['smiles'], new_rxn)

    return rxn

db = Database()
template = db.find('templates', {'name': 'temp1a'})
side_chains = db.find_all('side_chains')

for temp in template:
    for sc in side_chains:
        generate_rxn_temp(temp, sc, store=True)
