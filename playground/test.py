import json
from molecules import SideChainGenerator, DataInitializer, MonomerGenerator
from utils import create_logger, MongoDataBase
from time import time
from rdkit import Chem
from rdkit.Chem import Draw
# from logging import
# logger = create_logger(__name__)

# db = MongoDataBase()
# db.setup()

# init = DataInitializer()
# init.from_json()
# init.to_json()
# init.initialize_data()
# init.load_parent_side_chains()
# init.load_connections()
# init.load_rtemplates()
# init.load_monomers(group='modified_prolines', required=False)

# generator = SideChainGenerator()
# if generator.load_data():
#     start = time()
#     generator.generate()
#     print(time()-start)
#     generator.save_data()

generator = MonomerGenerator()
if generator.load_data():
    start = time()
    generator.generate_from_ids(['b0'], ['alpha'])
    print(time() - start)
    print(generator.result_data)

# generator = MonomerGenerator()
# if generator.load_data():
#     start = time()
#     generator.generate_parallel()
#     print(time() - start)
#     p_results = generator.result_data

# print(len(results))
# print(len(p_results))
