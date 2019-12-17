import exceptions
import os

from rdkit import Chem

import config
import importers
import initializers
import molecules
import project_io
import reactions


def connect_mols(*mols, map_nums, stereo=None, clear_map_nums=True):
    """
    Function for combining either one or two molecules at the positions marked by atom map numbers. This function also
    applies the specified stereochemistry at the connected position if applicable, and can clear the map numbers from
    the molecule after making the connection if desired.

    Args:
        map_nums (iterable): An iterable containing two integer elements that specify the map numbers on the atoms
            to connect.
        stereo (str, optional): The stereochemistry to place on the new connection. Can either be 'CW' or 'CCW'.
            Defaults to None.
        clear_map_nums (bool, optional): Whether to clear atom numbers after making the connection or not.
            Defaults to True.

    Raises:
        exceptions.MergeError: Raised if either no molecules or more than two are provided.

    Returns:
        RDKit Mol: The result of connecting the molecule(s) at the specified positions.
    """

    def update_hydrogen_counts(atom, clear_map_nums):
        """
        Inner method for clearing the atom map number and updating hydrogen counts.
        """

        if clear_map_nums:
            atom.SetAtomMapNum(0)

        if atom.GetSymbol() in ['N', 'O', 'S']:
            atom.SetNumExplicitHs(0)
        elif atom.GetSymbol() == 'C' and atom.GetNumExplicitHs() != 0:
            atom.SetNumExplicitHs(atom.GetTotalNumHs() - 1)

        return atom

    if len(mols) < 1 or len(mols) > 2:
        raise exceptions.MergeError('Can only merge 1 or 2 molecules at a time.')

    # find atoms that will form a bond together and update hydrogen counts
    combo, *mols = mols
    for mol in mols:
        combo = Chem.CombineMols(combo, mol)

    # find atoms that will form a bond together and update hydrogen counts
    combo = Chem.RWMol(combo)
    Chem.SanitizeMol(combo)
    try:
        atom1, atom2 = [update_hydrogen_counts(atom, clear_map_nums)
                        for atom in combo.GetAtoms() if atom.GetAtomMapNum() in map_nums]
    except ValueError:
        raise exceptions.MergeError('There must be exactly 2 map numbers across all molecules.')

    # create bond, remove hydrogens, and sanitize
    combo.AddBond(atom1.GetIdx(), atom2.GetIdx(), order=Chem.BondType.SINGLE)
    Chem.RemoveHs(combo)
    Chem.SanitizeMol(combo)

    # add stereochemistry as specified
    stereo_center = atom1 if atom1.GetHybridization() == Chem.HybridizationType.SP3 and atom1.GetTotalNumHs() != 2 else atom2
    if stereo == 'CCW':
        stereo_center.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
    elif stereo == 'CW':
        stereo_center.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)

    return Chem.MolFromSmiles(Chem.MolToSmiles(combo))


def create_regiosqm_smiles_file(sidechain_io):

    project_io.RawRegioSQMIO().save(filter(lambda x: x['connection'] == 'methyl', sidechain_io.load()))


def reset(data_format=config.DATA_FORMAT):

    if data_format == 'mongo':
        database = project_io.MongoDataBase()
        database.setup()

    initializers.RecordInitializer(data_format).initialize()
    importers.DataImporter(data_format).import_data()


def file_rotator(filepath):
    """
    Function that takes a filepath and attaches a underscore and a number before its extension to ensure uniqueness of
    the file name given by the filepath.

    Args:
        filepath (str): The path to the file.

    Returns:
        str: The augmented filepath.
    """

    idx = 0
    while True:
        new_fp = attach_file_num(filepath, idx)
        idx += 1
        if not (os.path.exists(new_fp) and os.path.isfile(new_fp)):
            return new_fp


def attach_file_num(filepath, file_num):
    """
    Function that inserts an underscore and the specified file number to the file name given in the filepath.

    Args:
        filepath (str): The filepath containing the file name to augment.
        file_num (int): The nile number to attach to the end of the file name in the filepath.

    Returns:
        str: The augmented filepath.
    """

    new_fp, ext = filepath.split('.')
    new_fp += '_' + str(file_num) + '.' + ext
    return new_fp


def get_file_num_range(filepath):
    """
    Function that scans the last directory in the given filepath for all files with the base file name specified in the
    filepath and determines the range of file numbers appended to the base file name that are already in use.

    Args:
        filepath (str): The path to the file.

    Returns:
        tuple[int]: A tuple of ints where the first argument is 0 and the second is the highest used file number of the
            base file name in the directory.
    """

    low = 0
    high = 0
    while True:
        new_fp = attach_file_num(filepath, high)
        if not (os.path.exists(new_fp) and os.path.isfile(new_fp)):
            return low, high
        high += 1


def atom_to_wildcard(atom):
    """
    Changes the atom to a wildcard atom.

    Args:
        atom (RDKit Atom): The atom to change to a wildcard.
    """

    atom.SetAtomicNum(0)
    atom.SetIsotope(0)
    atom.SetFormalCharge(0)
    atom.SetIsAromatic(False)
    atom.SetNumExplicitHs(0)


def get_templates():
    """
    Creates and returns a list of all template molecules.

    Returns:
        list: The template molecules
    """

    return [molecules.CinnamoylTemplate1(), molecules.CinnamoylTemplate2(), molecules.CinnamoylTemplate3()]


def get_connections():
    """
    Creates and returns a list of all connection molecules.

    Returns:
        list: The connection molecules.
    """

    return [molecules.EthylConnection()]


def get_backbones():
    """
    Creates and returns a list of all backbone moleucles.

    Returns:
        list: The backbone molecules.
    """

    return [molecules.AlphaBackBone(), molecules.Beta2BackBone(), molecules.Beta3BackBone()]


def get_hashed_molecules(func):
    """
    Closuer that returns a function that when called returns a dict of the molecules returned by func, where the keys
    are the names of the molecules and the values are the RDKit representation of the molecule.

    Args:
        func (function): A getter function that returns a list of Molecule objects such as get_templates().

    Returns:
        function: The function that will return the hashed molecules when called.
    """

    def hasher():
        hashed_molecules = {}
        for molecule in func():
            hashed_molecules[molecule.name] = molecule.mol
        return hashed_molecules

    return hasher


get_hashed_templates = get_hashed_molecules(get_templates)
get_hashed_connections = get_hashed_molecules(get_connections)
get_hashed_backbones = get_hashed_molecules(get_backbones)


def get_reactions():
    return [reactions.FriedelCrafts(), reactions.TsujiTrost(), reactions.PictetSpangler(),
            reactions.TemplatePictetSpangler(), reactions.PyrroloIndolene()]


def get_reactions_of_type(rxn_type):

    def reaction_getter():
        return filter(lambda x: rxn_type in x.type, get_reactions())

    return reaction_getter


get_bimolecular_reactions = get_reactions_of_type('bimolecular')
get_unimolecular_reactions = get_reactions_of_type('unimolecular')


def get_hashed_predictions(func):

    def hasher():
        hashed_predictions = {}
        for prediction in func():
            hashed_predictions[prediction['_id']] = prediction['predictions']
        return hashed_predictions

    return hasher


def get_regiosqm_predictions(data_format=config.DATA_FORMAT):

    if data_format == 'json':
        regio_io = project_io.JsonRegioSQMIO()
    elif data_format == 'mongo':
        regio_io = project_io.MongoRegioSQMIO()

    return regio_io.load()


get_hashed_regiosqm_predictions = get_hashed_predictions(get_regiosqm_predictions)
