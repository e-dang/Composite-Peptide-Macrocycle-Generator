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
            for sidechain in Chem.ReplaceSubstructs(parent_sidechain.mol, models.NON_ATTACHMENT_METHYL, connection.mol):
                new_sidechains.append(models.Sidechain.from_mol(sidechain, connection, parent_sidechain.shared_id))

        return new_sidechains
