from time import time

from rdkit import Chem

import config
import factories
import importers
import initializers
import planners
import project_io


def initialize():
    """
    Function for resetting the library by removing all data from the data/generated/ directory or the MongoDataBase
    depending on the data format specified in the config file, then importing the sidechains and monomers from the sdf
    files, initializing the id and index records, and then creating the SMILES files for RegioSQM and pKa predictions
    After the predictions have been completed, place them in the external directory within the data directory and call
    the import_predictions() method of the DataImporter class. This method should be called if predictions have not been
    made. Otherwise if you are looking to simply reset the library use the function reset().
    """

    if config.DATA_FORMAT == 'mongo':
        database = project_io.MongoDataBase()
        database.setup()

    initializers.RecordInitializer().initialize()
    importers.DataImporter().import_molecules()
    create_regiosqm_smiles_file()
    create_pka_smiles_file()


def reset():
    """
    Function for resetting the library by removing all data from the data/generated/ directory or the MongoDataBase
    depending on the data format specified in the config file, then importing the sidechains and monomers from the sdf
    files, as well as the predictions, and initializing the id and index records. This function assumes that the
    predictions are already in the data/external/ directory.
    """

    if config.DATA_FORMAT == 'mongo':
        database = project_io.MongoDataBase()
        database.setup()

    initializers.RecordInitializer().initialize()
    importers.DataImporter().import_data()


def create_regiosqm_smiles_file():
    """
    Function for creating the SMILES file needed as input by RegioSQM.
    """

    sidechain_io = project_io.get_sidechain_io()
    monomer_io = project_io.get_monomer_io()

    data = list(filter(lambda x: x['connection'] == 'methyl', sidechain_io.load()))
    data.extend(list(filter(lambda x: x['required'], monomer_io.load())))
    project_io.RawRegioSQMIO().save(data)


def create_pka_smiles_file():
    """
    Function for creating the SMILES file used for tabulating the pKa predictions.
    """

    data = []
    sidechain_io = project_io.get_sidechain_io()
    for sidechain in filter(lambda x: x['connection'] == 'methyl', sidechain_io.load()):
        atom_map = 1
        mol = Chem.Mol(sidechain['binary'])
        for atom in mol.GetAtoms():
            if atom.GetSymbol() in ('N', 'O', 'S') and atom.GetTotalNumHs() > 0:
                atom.SetAtomMapNum(atom_map)
                atom_map += 1

        if atom_map > 1:
            Chem.Kekulize(mol)
            data.append((sidechain['_id'], Chem.MolToSmiles(mol, kekuleSmiles=True)))

    project_io.RawpKaIO().save(data)


def run_serial(factory_arg, name):
    """
    Closure for creating a function that passes the given factory argument to a factory to be ran in a serial fashion
    and printing out associated data about the run such as total run time and number of data points generated. The
    function also returns the associated data about the run.

    Args:
        factory_arg (IFactoryArgument): An implementation of the IFactoryArgument class.
        name (str): The name of the data being generated, i.e. SideChains, Monomers, etc...

    Returns:
        func: A function that when called runs the factory with the given factory arguments serially.
    """

    def closure(**kwargs):
        factory = factories.MolFactory()
        start = time()
        factory.run_serial(factory_arg(**kwargs))
        print_str = f'{name}:', 'Time:', time() - start, 'Num Data:', factory.count
        print(*print_str)
        return print_str

    return closure


def run(factory_arg, name):
    """
    Closure for creating a function that passes the given factory argument to a factory to be ran in a parallel fashion
    and printing out associated data about the run such as total run time and number of data points generated. The
    function also returns the associated data about the run.

    Args:
        factory_arg (IFactoryArgument): An implementation of the IFactoryArgument class.
        name (str): The name of the data being generated, i.e. SideChains, Monomers, etc...

    Returns:
        func: A function that when called runs the factory with the given factory arguments serially.
    """

    def closure(**kwargs):
        factory = factories.MolFactory()
        start = time()
        factory.run(factory_arg(**kwargs))
        try:
            name_with_length = name + ' ' + str(kwargs['peptide_length'])
            print_str = f'{name_with_length}:', 'Time:', time() - start, 'Num Data:', factory.count
            print(*print_str)
            return print_str
        except KeyError:
            print(f'{name}:', 'Time:', time() - start, 'Num Data:', factory.count)

    return closure


run_sidechains = run_serial(factories.SCConnectionModificationArgs, 'Sidechains')
run_monomers = run_serial(factories.MonomerGenerationArgs, 'Monomers')
run_peptides = run(factories.PublicationPeptideGenerationArgs, 'Peptides')
run_template_peptides = run(factories.TemplatePeptideGenerationArgs, 'Template Peptides')
run_unimolecular_reactions = run(factories.UniMolecularReactionGenerationArgs, 'UniMolecular Reactions')
run_bimolecular_reactions = run(factories.BiMolecularReactionGenerationArgs, 'BiMolecular Reactions')
run_macrocycles = run(factories.MacrocycleGenerationArgs, 'Macrocycles')


def generate_peptide_plan(peptide_length, num_peptides):
    """
    Function for generating the peptide plan with a given number of peptides with the given peptide length. A peptide
    plan is a txt file with each line containing the indices of the monomers that will be used to generate the peptide.
    The is useful since it allows for tracking of which peptides have already been generated so that adding onto the
    library can be done without having to worry about duplicate peptides. The peptide plan is also generated to ensure
    each monomer is used at each possible position of a peptide at least once.

    Args:
        peptide_length (int): The desired length of the peptides.
        num_peptides (int): The number of peptides to generate. Must be at least as big as the total number of monomers
            * the peptide length.
    """

    monomer_io = project_io.get_monomer_io()
    planner = planners.PeptidePublicationPlanner(monomer_io, peptide_length, num_peptides)
    planner.create_plan()
