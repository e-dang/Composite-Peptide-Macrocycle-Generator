from rdkit import Chem
import argparse

def main(file, out):
    with open('/Users/ericdang/Documents/UCLA_Research/smiles/c_term/' + file, 'r') as fin, open('/Users/ericdang/Documents/UCLA_Research/smiles/amide_term/' + out, 'w') as fout:
        for line in fin:
            mol = Chem.MolFromSmiles(line)
            patt = Chem.MolFromSmarts('O[$(C(=O)C(*)NC(=O)*),$(C(=O)CNC(=O)*)]')    # pattern matches peptide backbone from C-term end

            # check if match was found
            if not mol.HasSubstructMatch(patt, useChirality=False):
                print('No match found for molecule: ' + line)   # something went wrong
                continue

            # match was found
            matches = mol.GetSubstructMatches(patt, useChirality=False)
            for atom_idx in matches[0]:     # do transform for only 1 match (possible multiple matches for aspartic acid)
                target = mol.GetAtomWithIdx(atom_idx)
                if  target.GetSymbol() == 'O':
                    Chem.Atom.SetAtomicNum(target, 7)   # change the oxygen to nitrogen

            Chem.SanitizeMol(mol)

            # write to out file
            fout.write(Chem.MolToSmiles(mol))
            fout.write('\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file', help='The text file containing smiles strings of template coupled peptides with c-term')
    parser.add_argument('out', help='The output text file')

    args = parser.parse_args()
    main(args.file, args.out)
