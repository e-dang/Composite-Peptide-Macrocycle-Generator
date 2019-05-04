from rdkit import Chem
from utils import read_mols
import argparse
from pathlib import Path

def main(fp_in, fp_out, sc_file, bb_file, out, stereo):
    side_chains = read_mols(sc_file, fp_in)
    backbone = read_mols(bb_file, fp_in)[0]

    # set atom map number for attachment point of backbone
    patt = Chem.MolFromSmarts('[NH2]C') # attachment point at alpha carbon
    matches = backbone.GetSubstructMatches(patt, useChirality=False)
    for pairs in matches:
        for atom_idx in pairs:
            atom = Chem.Mol.GetAtomWithIdx(backbone, atom_idx)
            if atom.GetSymbol() == 'C' and Chem.Atom.GetTotalNumHs(atom) == 2:
                atom.SetAtomMapNum(1)

    with open(fp_out + out, 'w') as f:
        for sc in side_chains:

            # set atom map number for attachment point of side chain
            patt = Chem.MolFromSmarts('C*') # attachment point at methyl carbon
            matches = sorted(sc.GetSubstructMatches(patt, useChirality=False), key=lambda x: x[0], reverse=True)

            if not len(matches):
                print("Couldn't find attchment point of:", Chem.MolToSmiles(sc))

            attach_complete = False # for checking if one atom pair per mol resulted in a successful attachment
            for pairs in matches:
                sc_atom_idx = None
                for atom_idx in pairs:
                    atom = Chem.Mol.GetAtomWithIdx(sc, atom_idx)
                    if atom.GetSymbol() == 'C' and Chem.Atom.GetTotalNumHs(atom) == 3 and not Chem.Atom.IsInRing(atom):
                        atom.SetAtomMapNum(2)
                        sc_atom_idx = atom_idx

                # combine side chain and backbone
                combo = Chem.RWMol(Chem.CombineMols(sc, backbone))

                # get reacting atom indicies
                bb_atom = None
                sc_atom = None
                for atom in combo.GetAtoms():
                    if atom.GetAtomMapNum() == 1:
                        bb_atom = atom.GetIdx()
                        atom.SetAtomMapNum(0)
                    elif atom.GetAtomMapNum() == 2:
                        sc_atom = atom.GetIdx()
                        atom.SetAtomMapNum(0)
                        Chem.Mol.GetAtomWithIdx(sc, sc_atom_idx).SetAtomMapNum(0)

                # create bond and adjust Stereochemistry
                try:
                    combo.AddBond(bb_atom, sc_atom, order=Chem.rdchem.BondType.SINGLE)

                    if stereo == 'S':
                        Chem.Mol.GetAtomWithIdx(combo, bb_atom).SetChiralTag(Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW)
                    else:
                        Chem.Mol.GetAtomWithIdx(combo, bb_atom).SetChiralTag(Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW)

                    Chem.SanitizeMol(combo)

                    # write smiles string to file
                    f.write(Chem.MolToSmiles(combo))
                    f.write('\n')

                    attach_complete = True
                except:
                    if attach_complete:
                        print('Could not create bond for side chain (success in previous atom pair):', Chem.MolToSmiles(sc))
                    else:
                        print('Could not create bond for side chain:', Chem.MolToSmiles(sc))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generates monomers by attaching side chains to a backbone and outputs result as a SMILES string. Sidechains and backbones are defined as .sdf files in ../chemdraw/ folder. Output text file is written to ../smiles/monomers/ folder')
    parser.add_argument('stereo', choices=['S', 'R'], help='Stereochemistry of the monomer.')
    parser.add_argument('-sc', dest='sc_file', default='sidechains_cust.sdf', help='The .sdf file containing monomer side chains.')
    parser.add_argument('-bb', dest='bb_file', default='monomer_backbone.sdf', help='The .sdf file containing monomer backbone.')
    parser.add_argument('-o', '--out', dest='out', default='custom_', help='The output text file.')
    parser.add_argument('-fi', '--fin', dest='fp_in', default='chemdraw', help='The input filepath relative to script')
    parser.add_argument('-fo', '--fout', dest='fp_out', default='smiles/monomers', help='The ouput filepath relative to script')

    args = parser.parse_args()

    # format out file argument based on Stereochemistry
    if args.out == 'custom_':
        args.out = args.out + args.stereo + '.txt'

    # get filepath to infile and outfile
    base_path = Path(__file__).resolve().parents[1] # filepath to base project directory
    fp_in = str(base_path / args.fp_in) + '/'
    fp_out = str(base_path / args.fp_out) + '/'

    main(fp_in, fp_out, args.sc_file, args.bb_file, args.out, args.stereo)
