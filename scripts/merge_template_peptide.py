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


def get_template(fp_in_t, temp_ind):
    """
    Retrieves the template based on the value of temp_ind passed into the function and passes the template to a
    helper function that modifies it in preparation for merging with peptide

    Args:
        fp_in_t (str): The filepath to the json file containing the template SMILES strings
        temp_ind (int): Index indicating which template to retrieve, 1 - temp1a, 2 - temp2, 3 - temp3

    Returns:
        rdkit Mol: The modified template
    """

    template = None
    with open(fp_in_t, 'r') as f:
        for ind, doc in enumerate(json.load(f)):
            if ind + 1 == temp_ind:
                smiles = doc['smiles']
                template = modify_template(Chem.MolFromSmiles(smiles), ind + 1)

    return template, smiles


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
    if ind == 1 or ind == 2:    # temp 1a or 1b
        patt1 = Chem.MolFromSmiles('O=C1CCC(=O)N1O')  # succinimide + extra oxygen atom
        patt2 = Chem.MolFromSmarts('[CH]=O')    # carbonyl that will form amide with peptide

    # remove reaction point substruct
    template_mod = AllChem.DeleteSubstructs(template, patt1)

    # find and set template peptide connection point
    matches = template_mod.GetSubstructMatches(patt2, useChirality=False)
    for pairs in matches:
        for atom_idx in pairs:
            atom = Chem.Mol.GetAtomWithIdx(template_mod, atom_idx)
            if atom.GetSymbol() == 'C' and Chem.Atom.GetTotalNumHs(atom) == 1:
                atom.SetAtomMapNum(1)

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

    # portion of peptide backbone containing n-term
    if n_term == ALPHA:
        patt = Chem.MolFromSmarts('[NH2]CC(=O)')
    elif n_term == BETA2 or n_term == BETA3:
        patt = Chem.MolFromSmarts('[NH2]CCC(=O)')

    # find n-term nitrogen and assign atom map number
    matches = peptide.GetSubstructMatches(patt, useChirality=False)
    for pairs in matches:
        for atom_idx in pairs:
            atom = Chem.Mol.GetAtomWithIdx(peptide, atom_idx)
            if atom.GetSymbol() == 'N' and Chem.Atom.GetTotalNumHs(atom) == 2:
                atom.SetAtomMapNum(2)

    # prep for modification
    combo = Chem.RWMol(Chem.CombineMols(peptide, template))

    # get reacting atom's indices in combo mol and remove atom map numbers
    pep_atom = None
    temp_atom = None
    for atom in combo.GetAtoms():
        if atom.GetAtomMapNum() == 1:
            temp_atom = atom.GetIdx()
            atom.SetAtomMapNum(0)
        elif atom.GetAtomMapNum() == 2:
            pep_atom = atom.GetIdx()
            atom.SetAtomMapNum(0)

    # create bond and sanitize
    combo.AddBond(temp_atom, pep_atom, order=Chem.rdchem.BondType.SINGLE)
    Chem.SanitizeMol(combo)

    return Chem.MolToSmiles(combo)


def main():
    parser = argparse.ArgumentParser(description='Connects each peptide defined in the input file to a template and '
                                     'write the resulting molecule to file as a SMILES string')
    parser.add_argument('-t', '--temp', dest='template', choices=[1, 2, 3], type=int, default=[1], nargs='+',
                        help='The template(s) to be used; 1 - temp1a,, 2 - temp2, 3 - temp3')
    parser.add_argument('-tin', '--temp_in', dest='temp_in', default='templates.json',
                        help='The json file containing template SMILES strings')
    parser.add_argument('-pin', '--pep_in', dest='pep_in', default=['length3_all.json'], nargs='+',
                        help='The json file(s) containing peptide SMILES strings')
    parser.add_argument('-o', '--out', dest='out', default=None, nargs='+',
                        help='The output json file(s) to write the resulting SMILES strings; default will out assign '
                        'file names')
    parser.add_argument('-fit', '--fin_t', dest='fp_in_t', default='smiles/templates',
                        help='The filepath to the template json file relative to the base project directory')
    parser.add_argument('-fip', '--fin_p', dest='fp_in_p', default='smiles/peptides',
                        help='The filepath to the peptide json file(s) relative to the base project directory')
    parser.add_argument('-fo', '--fout', dest='fp_out', default='smiles/template_peptide/c_term',
                        help='The filepath for the output json file relative to the base project directory')
    parser.add_argument('--show_progress', dest='progress', action='store_false',
                        help='Show progress bar. Defaults to False')

    args = parser.parse_args()

    # if output file(s) specified then check for proper number of specified files
    if args.out is not None:
        num_output = len(args.out)
        num_temp = len(args.template)
        num_pep = len(args.pep_in)
        if num_output != num_temp * num_pep:
            print('Number of output files needs to be equal to number of templates * number of peptide files')
            raise SystemExit

    # get full filepath to input files
    base_path = Path(__file__).resolve().parents[1]
    fp_in_t = str(base_path / args.fp_in_t / args.temp_in)
    fp_in_p = [str(base_path / args.fp_in_p / file) for file in args.pep_in]

    # connect all peptides in each peptide group to each template and write to correct output file
    template_names = ['_temp1', '_temp2', '_temp3']
    for i, peptide_fp in tqdm(enumerate(fp_in_p), disable=args.progress):
        for j, temp_ind in tqdm(enumerate(args.template), disable=args.progress):

            # create output filepath based on template and peptide lengths if output file(s) not specified
            outfile_name = args.pep_in[i].split('_')[0] + template_names[j] + '.json'
            fp_out = str(base_path / args.fp_out /
                         outfile_name) if args.out is None else str(args.out[i * len(args.template) + j])

            # prep corresponding template and peptides
            template, smiles = get_template(fp_in_t, temp_ind)

            # connect and write to file
            collection = []
            with open(fp_out, 'w') as fout, open(peptide_fp, 'r') as f_pep:
                for doc in tqdm(json.load(f_pep), disable=args.progress):
                    doc['template'] = smiles
                    doc['template-peptide'] = combine(template, Chem.MolFromSmiles(doc['peptide']), doc['N-term'])
                    del doc['N-term']
                    collection.append(doc)

                json.dump(collection, fout)
                print('Complete!')


if __name__ == '__main__':
    main()
