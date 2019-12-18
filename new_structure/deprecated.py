class SideChainGenerator(IGenerator):
    """
    Implementation of an IGenerator that takes a connection molecule and heterocycle and creates a set of
    sidechains by attaching the connection molecule to all valid positions on the heterocycle. Valid positions are
    determined in instace method is_valid_atom().
    """

    # atom map numbers
    CONNECTION_MAP_NUM = molecules.IConnectionMol.OLIGOMERIZATION_MAP_NUM
    HETEROCYCLE_MAP_NUM = 2
    MAP_NUMS = (CONNECTION_MAP_NUM, HETEROCYCLE_MAP_NUM)

    def get_args(self, data):

        return product(data.heterocycles, data.connections)

    def generate(self, args):

        self.sidechains = {}
        self.heterocycle, self.connection = args
        heterocycle = Chem.Mol(self.heterocycle['binary'])
        connection = self.connection.tagged_mol

        # check if connecting atom is atom mapped
        self.is_valid_connection(connection)

        # tag any methyl groups on heterocycle so they are not confused with methyl connection atoms in the monomer
        # generation step
        patt = Chem.MolFromSmarts('[CH3]')
        for i, atom in enumerate(chain.from_iterable(heterocycle.GetSubstructMatches(patt)), start=3):
            heterocycle.GetAtomWithIdx(atom).SetAtomMapNum(i)

        # make attachment at each atom
        for atom in heterocycle.GetAtoms():

            # detetmine atom eligibility
            if not self.is_valid_atom(atom):
                continue

            # merge parent side chain with conenction
            atom.SetAtomMapNum(self.HETEROCYCLE_MAP_NUM)
            sidechain = utils.connect_mols(heterocycle, connection, map_nums=self.MAP_NUMS)
            atom.SetAtomMapNum(0)

            # check for uniqueness and record results
            binary = sidechain.ToBinary()
            Chem.Kekulize(sidechain)
            self.sidechains[binary] = (Chem.MolToSmiles(sidechain, kekuleSmiles=True), atom.GetIdx())

        return self.format_data()

    def is_valid_atom(self, atom):
        """
        Helper method used to determine if a specific position on the heterocycle in valid.

        Args:
            atom (RDKit Atom): The candidate atom on the heterocycle that might have the connection molecule attached
                to it.

        Returns:
            bool: True if valid.
        """

        valid_carbon = atom.GetSymbol() == 'C' and 0 < atom.GetTotalNumHs() < 3
        valid_hetero = atom.GetSymbol() in ('N', 'O', 'S') and atom.GetTotalNumHs() != 0

        return valid_carbon or valid_hetero

    def is_valid_connection(self, connection):
        """
        Helper method that determines whether the connection molecule has an atom map number specifying the point of
        attachment or not.

        Args:
            connection (RDKit Mol): The connection molecule.

        Raises:
            exceptions.MissingMapNumber: Raised when there is no atom map number on the connection molecule.
        """

        if self.CONNECTION_MAP_NUM not in [atom.GetAtomMapNum() for atom in connection.GetAtoms()]:
            raise exceptions.MissingMapNumber('Connection molecule missing atom map number')

    def format_data(self):
        """
        Helper method that fills in a dict object for each new sidechain with the necessary data associated with that
        sidechain needed for record keeping.

        Returns:
            list: A list containing the newly created sidechain dicts.
        """

        data = []
        # chunk = len(sidechains) * self.connections.index(connection)
        chunk = len(self.sidechains) * 3
        for i, (binary, (kekule, conn_atom_idx)) in enumerate(self.sidechains.items()):
            data.append({'_id': self.heterocycle['_id'] + str(chunk + i),
                         'type': 'sidechain',
                         'binary': binary,
                         'kekule': kekule,
                         'conn_atom_idx': conn_atom_idx,
                         'heterocycle': self.heterocycle['_id'],
                         'connection': self.connection.name})

        return data


class JsonHeterocycleIO(AbstractJsonIO):
    """
    Implmentation of the AbstractJsonIO class for handling heterocycle data.
    """

    _FILEPATH = os.path.join(config.DATA_DIR, 'generated', 'heterocycles.json')

    def load(self):

        return super().from_json(self._FILEPATH)

    def save(self, data):

        super().to_json(self._FILEPATH, data)


class MongoHeterocycleIO(AbstractMongoIO):
    """
    Implmentation of the AbstractMongoIO class for handling heterocycle data.
    """

    _COLLECTION = config.COL1
    _QUERY = {'type': 'heterocycle'}

    def load(self):

        return super().from_mongo(self._COLLECTION, self._QUERY)

    def save(self, data):

        super().to_mongo(self._COLLECTION, data)
