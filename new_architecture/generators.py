import new_architecture.models as models
import new_architecture.utils as temp_utils
import macrocycles.utils as utils
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


class MonomerGenerator:
    MAP_NUMS = (models.Backbone.MAP_NUM, models.Sidechain.MAP_NUM)

    def generate(self, args):
        """
        Method that takes a sidechain and backbones creates new monomers.
        """

        sidechain, backbones = args
        sidechain_mol = sidechain.mapped_mol

        temp_utils.clear_isotopes(sidechain_mol)

        monomers = []
        for backbone in backbones:
            monomer = utils.connect_mols(sidechain_mol, backbone.mol, map_nums=self.MAP_NUMS)
            monomers.append(models.Monomer.from_mol(monomer, backbone, sidechain))

        return monomers
