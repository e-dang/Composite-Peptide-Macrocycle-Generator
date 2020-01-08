from time import time

from rdkit import Chem

import config
import factories
import importers
import initializers
import planners
import project_io


def initialize():

    if config.DATA_FORMAT == 'mongo':
        database = project_io.MongoDataBase()
        database.setup()

    initializers.RecordInitializer().initialize()
    importers.DataImporter().import_molecules()
    create_regiosqm_smiles_file()
    create_pka_smiles_file()


def reset():

    if config.DATA_FORMAT == 'mongo':
        database = project_io.MongoDataBase()
        database.setup()

    initializers.RecordInitializer().initialize()
    importers.DataImporter().import_data()


def create_regiosqm_smiles_file():

    sidechain_io = project_io.get_sidechain_io()
    monomer_io = project_io.get_monomer_io()

    data = list(filter(lambda x: x['connection'] == 'methyl', sidechain_io.load()))
    data.extend(list(filter(lambda x: x['required'], monomer_io.load())))
    project_io.RawRegioSQMIO().save(data)


def create_pka_smiles_file():

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
    def closure(**kwargs):
        factory = factory = factories.MolFactory()
        start = time()
        factory.run_serial(factory_arg(**kwargs))
        print(f'{name}:', 'Time:', time() - start, 'Num Data:', factory.count)

    return closure


def run(factory_arg, name):
    def closure(**kwargs):
        factory = factory = factories.MolFactory()
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

    monomer_io = project_io.get_monomer_io()
    planner = planners.PeptidePublicationPlanner(monomer_io, peptide_length, num_peptides)
    planner.create_plan()
