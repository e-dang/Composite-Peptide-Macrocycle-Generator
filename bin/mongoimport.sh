#!/bin/bash

mongoimport --db $1 --collection $2 --file ../smiles/$3 --jsonArray