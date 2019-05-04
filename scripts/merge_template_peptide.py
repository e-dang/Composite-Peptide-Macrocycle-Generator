from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

temp = Chem.SDMolSupplier('/Users/ericdang/Documents/UCLA_Research/chemdraw/templates.sdf')[0]

# remove reaction point substruct
patt = Chem.MolFromSmarts('O=C1CCC(=O)N1O')
temp = AllChem.DeleteSubstructs(temp, patt)


# find and set template peptide connection point
patt = Chem.MolFromSmarts('[CH]=O')
matches = temp.GetSubstructMatches(patt, useChirality=False)
for pairs in matches:
    for atom_idx in pairs:
        atom = Chem.Mol.GetAtomWithIdx(temp, atom_idx)
        if atom.GetSymbol() == 'C' and Chem.Atom.GetTotalNumHs(atom) == 1:
            atom.SetAtomMapNum(1)

# get peptides, find and set connection point
with open('/Users/ericdang/Documents/UCLA_Research/smiles/peptides/demo_10.txt', 'r') as fin, open('/Users/ericdang/Documents/UCLA_Research/smiles/template_peptide/demo_10.txt', 'w') as fout:
    patt = Chem.MolFromSmarts('NCC(=O)NCC(=O)')
    for line in fin.readlines():
        print(line)
        peptide = Chem.MolFromSmiles(line)
        matches = peptide.GetSubstructMatches(patt, useChirality=False)
        print(matches)
        for pairs in matches:
            for atom_idx in pairs:
                atom = Chem.Mol.GetAtomWithIdx(peptide, atom_idx)
                if atom.GetSymbol() == 'N' and Chem.Atom.GetTotalNumHs(atom) == 2:
                    atom.SetAtomMapNum(2)

        print(Chem.MolToSmiles(peptide))
        combo = Chem.RWMol(Chem.CombineMols(peptide, temp))

        pep_atom = None
        temp_atom = None
        for atom in combo.GetAtoms():
            if atom.GetAtomMapNum() == 1:
                temp_atom = atom.GetIdx()
                atom.SetAtomMapNum(0)
            elif atom.GetAtomMapNum() == 2:
                pep_atom = atom.GetIdx()
                atom.SetAtomMapNum(0)

        combo.AddBond(temp_atom, pep_atom, order=Chem.rdchem.BondType.SINGLE)
        Chem.SanitizeMol(combo)

        fout.write(Chem.MolToSmiles(combo))
        fout.write('\n')

        # Draw.MolToImage(combo).show()
