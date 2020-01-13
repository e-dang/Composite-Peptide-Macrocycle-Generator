import macrocycles.runners as runners

runners.reset()

runners.run_sidechains()
runners.run_monomers()
runners.run_unimolecular_reactions()
runners.run_bimolecular_reactions()
