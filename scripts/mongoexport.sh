#!/bin/bash

DATE=$(date +"%m%d%y")

mongoexport --db rxn_templates --collection candidates --out ~/Documents/UCLA_Research/json/${DATE}/candidates.json --pretty
mongoexport --db rxn_templates --collection reactions --out ~/Documents/UCLA_Research/json/${DATE}/reactions.json --pretty
mongoexport --db rxn_templates --collection side_chains --out ~/Documents/UCLA_Research/json/${DATE}/side_chains.json --pretty
mongoexport --db rxn_templates --collection templates --out ~/Documents/UCLA_Research/json/${DATE}/templates.json --pretty
