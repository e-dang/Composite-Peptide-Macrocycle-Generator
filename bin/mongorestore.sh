#!/bin/bash

mongoimport --db rxn_templates --collection reactions --file ../json/$1/rxn_templates/reactions.json --jsonArray --drop
mongoimport --db rxn_templates --collection regioSQM --file ../json/$1/rxn_templates/regioSQM.json --jsonArray --drop
mongoimport --db rxn_templates --collection side_chains --file ../json/$1/rxn_templates/side_chains.json --jsonArray --drop
mongoimport --db rxn_templates --collection templates --file ../json/$1/rxn_templates/templates.json --jsonArray --drop

mongoimport --db molecules --collection side_chains --file ../json/$1/molecules/side_chains.json --jsonArray --drop
mongoimport --db molecules --collection monomers --file ../json/$1/molecules/monomers.json --jsonArray --drop
mongoimport --db molecules --collection peptides --file ../json/$1/molecules/peptides.json --jsonArray --drop
mongoimport --db molecules --collection template_peptides --file ../json/$1/molecules/template_peptides.json --jsonArray --drop
mongoimport --db molecules --collection candidates --file ../json/$1/molecules/candidates.json --jsonArray --drop
mongoimport --db molecules --collection regio_filtered_candidates --file ../json/$1/molecules/regio_filtered_candidates.json --jsonArray --drop