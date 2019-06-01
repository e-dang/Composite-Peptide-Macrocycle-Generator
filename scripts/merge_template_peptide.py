import argparse
import json
from itertools import chain
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from tqdm import tqdm

ALPHA = 'alpha_amino_acid'
BETA2 = 'beta2_amino_acid'
BETA3 = 'beta3_amino_acid'

SUCCINIMIDE = 'O=C1CCC(=O)N1O'  # has extra oxygen
CARBONYL = '[CH]=O'

TEMP_ATOM_MAP_NUM = 1
PEP_ATOM_MAP_NUM = 2

TEMPLATE_NAMES = ['_temp1', '_temp2', '_temp3']


def get_templates(fp_temp, temp_inds):
    """
    Retrieves the template based on the value of temp_inds passed into the function and passes the template to a
    helper function that modifies it in preparation for merging with peptide

    Args:
        fp_temp (str): The filepath to the json file containing the template SMILES strings
        temp_inds (int): Indicies indicating which templates to retrieve; 1 - temp1, 2 - temp2, 3 - temp3

    Returns:
        list: Contains tuples, each of which contain the modified template as an rdkit Mol and the unmodified template
            SMILES string
    """

    template_data = []
    with open(fp_temp, 'r') as f:
        for temp_ind in temp_inds:
            for i, doc in enumerate(json.load(f), start=1):
                if i == temp_ind:
                    template_data.append((modify_template(Chem.MolFromSmiles(doc['smiles']), i), doc['smiles']))

    return template_data


def modify_template(template, ind):
    """
    Helper function of get_template(), which performs a deletion of the substructure corresponding to the leaving
    group upon amide bond formation and the assignment of an atom map number to reacting carbon atom

    Args:
        template (rdkit Mol): The rdkit Mol representation of the template
        ind (int): The index that idicates the identity of the template

    Returns:
        rdkit Mol: The modified template
    """

    # TODO: Need to implement same process for templates 2 and 3
    if ind == 1:    # temp 1
        patt1 = Chem.MolFromSmiles(SUCCINIMIDE)
        patt2 = Chem.MolFromSmarts(CARBONYL)    # carbonyl that will form amide with peptide

    # remove reaction point substruct
    template_mod = AllChem.DeleteSubstructs(template, patt1)

    # find and set template peptide connection point
    matches = template_mod.GetSubstructMatches(patt2, useChirality=False)
    for pairs in matches:
        for atom_idx in pairs:
            atom = Chem.Mol.GetAtomWithIdx(template_mod, atom_idx)
            if atom.GetSymbol() == 'C' and Chem.Atom.GetTotalNumHs(atom) == 1:
                atom.SetAtomMapNum(TEMP_ATOM_MAP_NUM)

    return template_mod


def combine(template, peptide, n_term):
    """
    Convert the peptide SMILES string to an rdkit Mol and combine it with the template through an amide linkage and the
    designated reacting site

    Args:
        template (rdkit Mol): The template with pre-set atom map number for reacting atom
        peptide (rdkit Mol): The peptide to be attached to the template

    Returns:
        str: The SMILES string of the molecule resulting from merging the template with the peptide
    """

    # print('template', template)
    # print('peptide', peptide)
    # # portion of peptide backbone containing n-term
    # if n_term == ALPHA:
    #     patt = Chem.MolFromSmarts('[NH2]CC(=O)')
    # elif n_term in (BETA2, BETA3):
    #     patt = Chem.MolFromSmarts('[NH2]CCC(=O)')

    # find n-term nitrogen and assign atom map number
    matches = peptide.GetSubstructMatches(Chem.MolFromSmarts('[NH2]'), useChirality=False)
    # print('matches', matches)
    for pairs in matches:
        for atom_idx in pairs:
            atom = Chem.Mol.GetAtomWithIdx(peptide, atom_idx)
            if atom.GetSymbol() == 'N' and Chem.Atom.GetTotalNumHs(atom) == 2:  # primary amines only
                atom.SetAtomMapNum(PEP_ATOM_MAP_NUM)

    # prep for modification
    combo = Chem.RWMol(Chem.CombineMols(peptide, template))

    # get reacting atom's indices in combo mol and remove atom map numbers
    pep_atom = None
    temp_atom = None
    for atom in combo.GetAtoms():
        if atom.GetAtomMapNum() == TEMP_ATOM_MAP_NUM:
            temp_atom = atom.GetIdx()
            atom.SetAtomMapNum(0)
        elif atom.GetAtomMapNum() == PEP_ATOM_MAP_NUM:
            pep_atom = atom.GetIdx()
            atom.SetAtomMapNum(0)

    # print('pep_atom', pep_atom)
    # print('temp_atom', temp_atom)
    # create bond
    combo.AddBond(temp_atom, pep_atom, order=Chem.rdchem.BondType.SINGLE)

    try:
        Chem.SanitizeMol(combo)
    except ValueError:
        print('Error!')
        print('Template:', Chem.MolToSmiles(template))
        print('Peptide:', Chem.MolToSmiles(peptide) + '\n')

    return Chem.MolToSmiles(combo)


def merge_templates_peptides(temp_inds, fp_temp, fp_peps, fp_out, show_progress=False):

    template_data = get_templates(fp_temp, temp_inds)
    for i, template in tqdm(enumerate(template_data), desc='Templates:', disable=show_progress):
        for fp_pep in tqdm(fp_peps, desc='Peptide Files:', disable=show_progress):

            # create output filepath based on template and peptide lengths
            outfile = fp_pep.split('/')[-1].strip('.json') + TEMPLATE_NAMES[i] + '.json'

            merged_data = merge_template_peptide(template, fp_pep, show_progress)
            write_data(merged_data, fp_out, outfile)


def merge_template_peptide(template, fp_pep, show_progress=False):

    collection = []
    with open(fp_pep, 'r') as f:
        for doc in tqdm(json.load(f), desc='Peptides:', disable=show_progress):
            doc['template'] = template[1]
            doc['template-peptide'] = combine(template[0], Chem.MolFromSmiles(doc['peptide']), doc['N-term'])
            del doc['N-term']
            collection.append(doc)

    return collection


def write_data(merged_data, fp_out, outfile):

    fp_out = str(fp_out / outfile)

    with open(fp_out, 'w') as f:
        json.dump(merged_data, f)

    print('Complete!')
    return True


def main():
    parser = argparse.ArgumentParser(description='Connects each peptide defined in the input file to a template and '
                                     'write the resulting molecule to file as a SMILES string')
    parser.add_argument('templates', choices=[1, 2, 3], type=int, nargs='+',
                        help='The template(s) to be used; 1 - temp1,, 2 - temp2, 3 - temp3')
    parser.add_argument('--f_temp', dest='f_temp', default='templates.json',
                        help='The json file containing template SMILES strings')
    parser.add_argument('--f_pep', dest='f_pep', default=['length3_all.json'], nargs='+',
                        help='The json file(s) containing peptide SMILES strings')
    parser.add_argument('--fp_temp', dest='fp_temp', default='smiles/templates',
                        help='The filepath to the template json file relative to the base project directory')
    parser.add_argument('--fp_pep', dest='fp_pep', default='smiles/peptides',
                        help='The filepath to the peptide json file(s) relative to the base project directory')
    parser.add_argument('--fp_out', dest='fp_out', default='smiles/template_peptide/c_term',
                        help='The filepath for the output json file relative to the base project directory')
    parser.add_argument('--show_progress', dest='progress', action='store_false',
                        help='Show progress bar. Defaults to False')

    args = parser.parse_args()

    # get full filepath to input files
    base_path = Path(__file__).resolve().parents[1]
    fp_temp = str(base_path / args.fp_temp / args.f_temp)
    fp_peps = [str(base_path / args.fp_pep / file) for file in args.f_pep]
    fp_out = base_path / args.fp_out

    merge_templates_peptides(args.templates, fp_temp, fp_peps, fp_out, args.progress)


if __name__ == '__main__':
    main()
