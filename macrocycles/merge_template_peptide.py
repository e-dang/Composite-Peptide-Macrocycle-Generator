import argparse
import json
from copy import deepcopy
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem
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
    # if ind == 1:    # temp 1
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
        list: Contains the SMILES strings of the molecule resulting from merging the template with the peptide
    """

    # portion of peptide backbone containing n-term
    if n_term == ALPHA:
        patt = Chem.MolFromSmarts('[NH;R]CC(=O)')
    elif n_term in (BETA2, BETA3):
        patt = Chem.MolFromSmarts('[NH;R]CCC(=O)')

    # find any primary amine (not amide) or proline n-terminus
    matches = list(peptide.GetSubstructMatches(Chem.MolFromSmarts('[$([NH2]);!$([NH2]C(=O)*)]'), useChirality=False))
    matches.extend(list(peptide.GetSubstructMatches(patt, useChirality=False)))

    # for each eligible nitrogen form a connection
    results = []
    for pairs in matches:

        # assign atom map number
        mapped_atom = None
        for atom_idx in pairs:
            atom = Chem.Mol.GetAtomWithIdx(peptide, atom_idx)
            if atom.GetSymbol() == 'N':  # primary amines only
                atom.SetAtomMapNum(PEP_ATOM_MAP_NUM)
                mapped_atom = atom
                break

        # prep for modification
        combo = Chem.RWMol(Chem.CombineMols(deepcopy(peptide), deepcopy(template)))
        mapped_atom.SetAtomMapNum(0)

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

        # merge
        try:
            combo.AddBond(temp_atom, pep_atom, order=Chem.rdchem.BondType.SINGLE)
            Chem.SanitizeMol(combo)
            results.append(Chem.MolToSmiles(combo))
        except ValueError:
            print('Sanitize Error!')
            print('Template:', Chem.MolToSmiles(template))
            print('Peptide:', Chem.MolToSmiles(peptide) + '\n')
        except Exception:
            print('Bond Addition Error!')
            print('Template:', Chem.MolToSmiles(template))
            print('Peptide:', Chem.MolToSmiles(peptide) + '\n')

    return results


def merge_templates_peptides(temp_inds, fp_temp, fp_peps, show_progress=False):
    """
    Top level function that calls merge_template_peptide() with every combination of specified templates and peptide
    files.

    Args:
        temp_inds (list): Contains ints which specify which templates to use
        fp_temp (str): The filepath to the file containing the templates
        fp_peps (list): Contains the filepath(s) to the file(s) containing the peptides
        show_progress (bool, optional): Show progress bar. Defaults to False.

    Yields:
        list: A list containing dictionaries of the data associated with the merged template_peptide molecule
        str: The corresponding output file name for the data
    """

    template_data = get_templates(fp_temp, temp_inds)
    for temp_ind, template in tqdm(zip(temp_inds, template_data), desc='Templates:', disable=show_progress):
        for fp_pep in tqdm(fp_peps, desc='Peptide Files:', disable=show_progress):

            # create output filepath based on template and peptide lengths
            outfile = fp_pep.split('/')[-1].strip('.json') + TEMPLATE_NAMES[temp_ind - 1] + '.json'

            template_peptide = merge_template_peptide(template, fp_pep, show_progress)
            yield template_peptide, outfile


def merge_template_peptide(template, fp_pep, show_progress=False):
    """
    Helper function to merge_templates_peptides(). Gets peptide data from file and calls combine() on each peptide -
    template combination.

    Args:
        template (tuple): Contains the template as an rdkit Mol and the SMILES string of the unmodified template
        fp_pep (str): The filepath to the file containing the peptides
        show_progress (bool, optional): Show progress bar. Defaults to False.

    Returns:
        list: A list containing dictionaries of the associated merged template_peptide data
    """

    collection = []
    with open(fp_pep, 'r') as f:
        for doc in tqdm(json.load(f), desc='Peptides:', disable=show_progress):
            doc['template'] = template[1]
            doc['template-peptide'] = combine(template[0], Chem.MolFromSmiles(doc['peptide']), doc['N-term'])
            del doc['N-term']
            collection.append(doc)

    return collection


def write_data(merged_data, fp_out):
    """
    Writes the merged data to a json file.

    Args:
        merged_data (iterable): An containing dictionaries of the associated merged template_peptide data and the
            corresponding output file name
        fp_out (str): The incomplete filepath to the output file. Filepath is completed inside this function

    Returns:
        bool: True if successful
    """
    for template_peptide, outfile in merged_data:
        fp_out = str(fp_out / outfile)

        with open(fp_out, 'w') as f:
            json.dump(template_peptide, f)

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

    merged_data = merge_templates_peptides(args.templates, fp_temp, fp_peps, args.progress)
    write_data(merged_data, fp_out)


if __name__ == '__main__':
    main()
