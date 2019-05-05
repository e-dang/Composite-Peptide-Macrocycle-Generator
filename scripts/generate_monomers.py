import argparse
from pathlib import Path
from tqdm import tqdm

from rdkit import Chem

from utils import read_mols


def get_backbone(fp_in_bb):
    """
    Gets the designated backbone structure from file and designates the attachment point of the
    side chain with an atom map number

    Args:
        fp_in_bb (string): The filepath to the file containing monomer backbone structures

    Returns:
        rdkit Mol: The backbone molecule
    """

    backbone = read_mols(fp_in_bb)[0]

    # TODO: allow for creation of monomers using different backbones
    # set atom map number for attachment point of backbone
    patt = Chem.MolFromSmarts('[NH2]C')  # attachment point at alpha carbon
    matches = backbone.GetSubstructMatches(patt, useChirality=False)
    for pairs in matches:
        for atom_idx in pairs:
            atom = Chem.Mol.GetAtomWithIdx(backbone, atom_idx)
            if atom.GetSymbol() == 'C' and Chem.Atom.GetTotalNumHs(atom) == 2:
                atom.SetAtomMapNum(1)

    return backbone


def get_side_chains(fp_in_sc):
    """
    Gets all side chain structures from all files and loads them into an array.

    Args:
        fp_in_sc (string): The filepath to the file(s) containing the desired side chains

    Returns:
        list: A list containing the side chains as rdkit Mols
    """

    side_chains = []
    for file in fp_in_sc:
        side_chains.extend([mol for mol in read_mols(file)])

    return side_chains


def create_monomer(backbone, side_chain, stereo):
    """
    Connects the side chain to the backbone at the designated attachment points.

    Args:
        backbone (rdkit Mol): The backbone to which the side chain will be attached to
        side_chain (rdkit Mol): The side chain being attached to the backbone structure

    Returns:
        rdkit Mol: The resulting monomer from connecting the backbone and side chain molecules
    """

    # set atom map number for attachment point of side chain
    # TODO: allow for attachment to longer alkyl chains such as ethyl carbons
    patt = Chem.MolFromSmarts('C*')  # attachment point at methyl carbon
    matches = sorted(side_chain.GetSubstructMatches(patt, useChirality=False), key=lambda x: x[0], reverse=True)

    if not len(matches):
        print("Couldn't find attchment point of:", Chem.MolToSmiles(side_chain))

    attach_complete = False  # for checking if one atom pair per mol resulted in a successful attachment
    monomers = []
    for pairs in matches:
        sc_atom_idx = None
        for atom_idx in pairs:
            atom = Chem.Mol.GetAtomWithIdx(side_chain, atom_idx)
            if atom.GetSymbol() == 'C' and Chem.Atom.GetTotalNumHs(atom) == 3 and not Chem.Atom.IsInRing(atom):
                atom.SetAtomMapNum(2)
                sc_atom_idx = atom_idx

        # combine side chain and backbone
        combo = Chem.RWMol(Chem.CombineMols(side_chain, backbone))

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
                Chem.Mol.GetAtomWithIdx(side_chain, sc_atom_idx).SetAtomMapNum(0)

        # create bond and adjust Stereochemistry
        try:
            combo.AddBond(bb_atom, sc_atom, order=Chem.rdchem.BondType.SINGLE)

            if stereo == 'S':
                Chem.Mol.GetAtomWithIdx(combo, bb_atom).SetChiralTag(Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW)
            else:
                Chem.Mol.GetAtomWithIdx(combo, bb_atom).SetChiralTag(Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW)

            Chem.SanitizeMol(combo)
            monomers.append(combo)
            attach_complete = True
        except:
            if attach_complete:
                print('Could not create bond for side chain (success in previous atom pair):',
                      Chem.MolToSmiles(side_chain))
            else:
                print('Could not create bond for side chain:', Chem.MolToSmiles(side_chain))

    return monomers


def main():
    """
    Parses arguments and passes them to functions to create monomers and write the resulting SMILES strings to file.
    """

    parser = argparse.ArgumentParser(
        description='Generates monomers by attaching side chains to a backbone and outputs result as a SMILES string. '
        'Sidechains and backbones are defined as .sdf files in ../chemdraw/ folder. Output text file is '
        'written to ../smiles/monomers/ folder')
    parser.add_argument('stereo', choices=['S', 'R'], help='Stereochemistry of the monomer.')
    parser.add_argument('-sc', dest='sc_file', nargs='+',
                        default=['sidechains_cust.sdf'], help='The sdf file containing monomer side chains.')
    parser.add_argument('-bb', dest='bb_file', default='monomer_backbone.sdf',
                        help='The .sdf file containing monomer backbone.')
    parser.add_argument('-o', '--out', dest='out', default='custom_', help='The output text file.')
    parser.add_argument('-fi', '--fin', dest='fp_in', default='chemdraw', help='The input filepath relative to script')
    parser.add_argument('-fo', '--fout', dest='fp_out', default='smiles/monomers',
                        help='The ouput filepath relative to script')

    args = parser.parse_args()

    # format outfile argument based on Stereochemistry
    if args.out == 'custom_':
        args.out = args.out + args.stereo + '.txt'

    # get filepath to infiles and outfile
    base_path = Path(__file__).resolve().parents[1]  # filepath to base project directory
    fp_in_sc = [str(base_path / args.fp_in / file) for file in args.sc_file]
    fp_in_bb = str(base_path / args.fp_in / args.bb_file)
    fp_out = str(base_path / args.fp_out / args.out)

    # get backbone and all side chains
    backbone = get_backbone(fp_in_bb)
    side_chains = get_side_chains(fp_in_sc)

    # create monomer and write to out file
    with open(fp_out, 'w') as f:
        for side_chain in tqdm(side_chains):
            monomers = create_monomer(backbone, side_chain, args.stereo)

            for monomer in monomers:
                f.write(Chem.MolToSmiles(monomer))
                f.write('\n')


if __name__ == '__main__':
    main()
