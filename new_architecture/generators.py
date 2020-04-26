import new_architecture.models as models
from rdkit import Chem


class SidechainModifier:

    def generate(self, args):
        """
        Method for creating new sidechain molecules from other sidechains by changing the methyl connection to a
        different connection type.
        """

        new_sidechains = []
        parent_sidechain, connections = args

        # replace the designated attachment position with each type of connection
        for connection in connections:
            for sidechain in Chem.ReplaceSubstructs(parent_sidechain.mol, models.SC_ATTACHMENT_POINT, connection.mol):
                new_sidechains.append(models.Sidechain.from_mol(sidechain, connection, parent_sidechain.shared_id))

        return new_sidechains


# class MonomerGenerator:
#     BB_MAP_NUM = molecules.IBackBoneMol.OLIGOMERIZATION_MAP_NUM
#     SC_MAP_NUM = 2
#     MAP_NUMS = (BB_MAP_NUM, SC_MAP_NUM)

#     def generate(self, args):
#         """
#         Method that takes a sidechain and uses it to create new monomer molecules by attaching it to each type of
#         backbone molecule listed in molecules.py.
#         """

#         monomers = []
#         sidechain, backbones = args
#         sidechain_mol = sidechain.mol

#         # find candidate attachment point at terminal end of alkyl chains
#         matches = sidechain.GetSubstructMatches(Chem.MolFromSmarts('[CH3;!13CH3]'))
#         if len(matches) != 1:
#             raise exceptions.InvalidMolecule(
#                 f'Sidechain {Chem.MolToSmiles(sidechain)} is missing an attachment point')

#         # set atom map number for attachment point of sidechain
#         for atom_idx in chain.from_iterable(matches):
#             sidechain.GetAtomWithIdx(atom_idx).SetAtomMapNum(self.SC_MAP_NUM)

#         utils.clear_isotopes(sidechain)

#         # connect monomer and backbone
#         for backbone in backbones:
#             monomer = utils.connect_mols(sidechain.mol, backbone.mol, map_nums=self.MAP_NUMS)

#         return self.format_data()
