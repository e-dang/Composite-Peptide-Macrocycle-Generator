from rdkit import Chem
from utils import read_mols, Database
from rdkit.Chem import Draw

side_chains = read_mols('sidechains_cust.sdf')
db = Database()
for sc in side_chains:

    # print(Chem.MolToSmiles(sc))
    # print('\n')
    smiles = Chem.MolToSmiles(sc)
    # find atom map number for attachment point of side chain
    patt = Chem.MolFromSmarts('[CH3]*')
    matches = sc.GetSubstructMatches(patt, useChirality=False)
    for pairs in matches:

        old_atom_idx = None
        for atom_idx in pairs:
            atom = Chem.Mol.GetAtomWithIdx(sc, atom_idx)
            if atom.GetSymbol() == 'C' and Chem.Atom.GetTotalNumHs(atom) == 3:
                atom.SetAtomMapNum(4)   # reaction template side chain attachment point index
                old_atom_idx = atom_idx

        # set atom map numbers for every possible reacting atom
        for atom in sc.GetAtoms():
            valid = False
            if atom.GetSymbol() == 'C' and atom.GetAtomMapNum() != 4 and Chem.Atom.GetTotalNumHs(atom) != 0 and Chem.Atom.GetIsAromatic(atom):
                atom.SetAtomMapNum(2)   # reaction template side chain reacting index
                valid = True
            elif atom.GetSymbol() == 'N' and Chem.Atom.GetTotalNumHs(atom) != 0:
                atom.SetAtomMapNum(2)
                valid = True
            elif atom.GetSymbol() == 'O' and Chem.Atom.GetTotalNumHs(atom) == 1:
                atom.SetAtomMapNum(2)
                valid = True

            if valid:
                mol_str = Chem.MolToSmiles(sc)
                ind = mol_str.find('[CH3:4]') + 7 # length of atom map string
                mol_str = mol_str[:ind] + '[*:4]' + mol_str[ind:]
                mol_str = mol_str.replace('[CH3:4]', 'C')

                atom.SetAtomMapNum(0)

                # print(mol_str)
                # Draw.MolToImage(Chem.MolFromSmarts(mol_str)).show()

                try:
                    Chem.MolFromSmarts(mol_str)
                except:
                    print('Cant convert smarts string to mol:')
                    print('smiles: ' + smiles + '\nsmarts: ' + mol_str)

                db.insert_sidechain(mol_str, smiles, 4, 2)

        Chem.Mol.GetAtomWithIdx(sc, old_atom_idx).SetAtomMapNum(0)
