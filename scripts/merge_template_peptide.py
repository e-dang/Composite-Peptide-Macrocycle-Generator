import argparse
from itertools import chain
from pathlib import Path
from tqdm import tqdm

from rdkit import Chem
from rdkit.Chem import AllChem, Draw


def get_template(fp_in_t, temp_ind):

    template = None
    with open(fp_in_t, 'r') as f:
        for ind, smiles in enumerate(f.readlines()):
            if ind + 1 == temp_ind:
                template = modify_template(Chem.MolFromSmiles(smiles), ind + 1)

    return template


def modify_template(template, ind):

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


def get_peptides(peptide_fp):

    with open(peptide_fp, 'r') as f:
        return f.readlines()


def combine(template, peptide):

    # portion of peptide backbone containing n-term
    patt = Chem.MolFromSmarts('NCC(=O)NCC(=O)')

    # find n-term nitrogen and assign atom map number
    peptide = Chem.MolFromSmiles(peptide)
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
    parser.add_argument('-t', '--temp', dest='template', choices=[1, 2, 3, 4], type=int, default=[1], nargs='+',
                        help='The template(s) to be used')
    parser.add_argument('-tin', '--temp_in', dest='temp_in', default='templates.txt',
                        help='The text file containing template SMILES strings')
    parser.add_argument('-pin', '--pep_in', dest='pep_in', default=['length3_all.txt'], nargs='+',
                        help='The text file(s) containing peptide SMILES strings')
    parser.add_argument('-o', '--out', dest='out', default=None, nargs='+',
                        help='The output text file(s) to write the resulting SMILES strings; default will out assign '
                        'file names')
    parser.add_argument('-fit', '--fin_t', dest='fp_in_t', default='smiles/templates',
                        help='The filepath to the template text file relative to the base project directory')
    parser.add_argument('-fip', '--fin_p', dest='fp_in_p', default='smiles/peptides',
                        help='The filepath to the peptide text file(s) relative to the base project directory')
    parser.add_argument('-fo', '--fout', dest='fp_out', default='smiles/template_peptide/c_term',
                        help='The filepath for the output text file relative to the base project directory')

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
    template_names = ['_temp1a', '_temp1b', '_temp2', '_temp3']
    for i, peptide_fp in tqdm(enumerate(fp_in_p)):
        for j, temp_ind in tqdm(enumerate(args.template)):

            # create output filepath based on template and peptide lengths if output file(s) not specified
            outfile_name = args.pep_in[i].split('_')[0] + template_names[j] + '.txt'
            fp_out = str(base_path / args.fp_out /
                         outfile_name) if args.out is None else str(args.out[i * len(args.template) + j])

            # prep corresponding template and peptides
            template = get_template(fp_in_t, temp_ind)
            peptides = get_peptides(peptide_fp)

            # connect and write to file
            with open(fp_out, 'w') as f:
                for peptide in tqdm(peptides):
                    f.write(combine(template, peptide))
                    f.write('\n')


if __name__ == '__main__':
    main()
