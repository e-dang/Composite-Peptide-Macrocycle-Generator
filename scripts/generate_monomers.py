import argparse
import json
from pathlib import Path

from rdkit import Chem
from tqdm import tqdm

from utils import read_mols


def get_backbone(fp_in_bb):
    """
    Gets all backbone structures from file

    Args:
        fp_in_bb (string): The filepath to the file containing monomer backbone structures

    Returns:
        list: A list of tuples containing the backbones as (rdkit Mol, type of backbone)
    """

    backbones = []
    with open(fp_in_bb, 'r') as f:
        for doc in json.load(f):
            backbones.append((Chem.MolFromSmarts(doc['backbone']), doc['type']))

    return backbones


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
        with open(file, 'r') as f:
            for doc in json.load(f):
                side_chains.append(Chem.MolFromSmiles(doc['side_chain']))

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
    patt = Chem.MolFromSmarts('[CH3]')  # attachment point at methyl carbon
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

            if stereo == 'CCW':
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
        'Need to specify stereochemistry and whether the set of monomers being generated should be added to the set of '
        'required monomers for a peptide. Sidechains and backbones are defined as json files in ../smiles/pre_monomer/ '
        'folder. Output json file is written to ../smiles/monomers/ folder')
    parser.add_argument('stereo', choices=['CW', 'CCW'], help='Stereochemistry of the monomer.')
    parser.add_argument('-sc', dest='sc_file', nargs='+',
                        default=['side_chains.json'], help='The json file containing monomer side chains.')
    parser.add_argument('-bb', dest='bb_file', default='backbones.json',
                        help='The json file containing monomer backbone.')
    parser.add_argument('-o', '--out', dest='out', default='custom_', help='The output json file.')
    parser.add_argument('-fi', '--fin', dest='fp_in', default='smiles/pre_monomer',
                        help='The input filepath relative to script')
    parser.add_argument('-fo', '--fout', dest='fp_out', default='smiles/monomers',
                        help='The ouput filepath relative to script')
    parser.add_argument('--show_progress', dest='progress', action='store_false',
                        help='Show progress bar. Defaults to False')
    parser.add_argument('--required', dest='required', action='store_true',
                        help='Designates whether these sets of monomers should be added to the set of required '
                        'monomers for a peptide. Defaults to False')

    args = parser.parse_args()

    # format outfile argument based on Stereochemistry
    if args.out == 'custom_':
        args.out = args.out + args.stereo + '.json'

    # get filepath to infiles and outfile
    base_path = Path(__file__).resolve().parents[1]  # filepath to base project directory
    fp_in_sc = [str(base_path / args.fp_in / file) for file in args.sc_file]
    fp_in_bb = str(base_path / args.fp_in / args.bb_file)
    fp_out = str(base_path / args.fp_out / args.out)

    # get backbone and all side chains
    backbones = get_backbone(fp_in_bb)
    side_chains = get_side_chains(fp_in_sc)

    # create monomer
    collection = []
    for backbone in tqdm(backbones, disable=args.progress):
        for side_chain in tqdm(side_chains, disable=args.progress):
            monomers = create_monomer(backbone[0], side_chain, args.stereo)

            for monomer in monomers:
                doc = {}
                doc['monomer'] = Chem.MolToSmiles(monomer)
                doc['type'] = backbone[1]
                doc['stereo'] = args.stereo
                doc['side_chain'] = Chem.MolToSmiles(side_chain)
                doc['required'] = args.required
                collection.append(doc)

    # write to file
    with open(fp_out, 'w') as f:
        json.dump(collection, f)


if __name__ == '__main__':
    main()
