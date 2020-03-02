import macrocycles.runners as runners
from time import sleep

runners.reset()

runners.run_sidechains()
runners.run_monomers()
sleep(1)
runners.run_unimolecular_reactions()
sleep(1)
runners.run_bimolecular_reactions()
