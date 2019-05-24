#!/bin/bash

DATE=$(date +"%m%d%y")

mongoexport --db rxn_templates --collection reactions --out ~/Documents/UCLA_Research/macrocycles/json/${DATE}/rxn_templates/reactions.json --pretty --jsonArray
mongoexport --db rxn_templates --collection side_chains --out ~/Documents/UCLA_Research/macrocycles/json/${DATE}/rxn_templates/side_chains.json --pretty --jsonArray
mongoexport --db rxn_templates --collection templates --out ~/Documents/UCLA_Research/macrocycles/json/${DATE}/rxn_templates/templates.json --pretty --jsonArray
mongoexport --db rxn_templates --collection regioSQM --out ~/Documents/UCLA_Research/macrocycles/json/${DATE}/rxn_templates/regioSQM.json --pretty --jsonArray

mongoexport --db molecules --collection candidates --out ~/Documents/UCLA_Research/macrocycles/json/${DATE}/molecules/candidates.json --pretty --jsonArray
mongoexport --db molecules --collection monomers --out ~/Documents/UCLA_Research/macrocycles/json/${DATE}/molecules/monomers.json --pretty --jsonArray
mongoexport --db molecules --collection peptides --out ~/Documents/UCLA_Research/macrocycles/json/${DATE}/molecules/peptides.json --pretty --jsonArray
mongoexport --db molecules --collection side_chains --out ~/Documents/UCLA_Research/macrocycles/json/${DATE}/molecules/side_chains.json --pretty --jsonArray
mongoexport --db molecules --collection template_peptides --out ~/Documents/UCLA_Research/macrocycles/json/${DATE}/molecules/template_peptides.json --pretty --jsonArray
mongoexport --db molecules --collection regio_filtered_candidates --out ~/Documents/UCLA_Research/macrocycles/json/${DATE}/molecules/regio_filtered_candidates.json --pretty --jsonArray