import runners

NUM_PEPTIDES = 200000
PEPTIDE_LENGTH_LOW = 3
PEPTIDE_LENGTH_HIGH = 5

# create peptide plan
for i in range(PEPTIDE_LENGTH_LOW, PEPTIDE_LENGTH_HIGH + 1):
    runners.generate_peptide_plan(i, NUM_PEPTIDES)

# create peptides
for i in range(PEPTIDE_LENGTH_LOW, PEPTIDE_LENGTH_HIGH + 1):
    runners.run_peptides(peptide_length=i)

# create template_peptides
for i in range(PEPTIDE_LENGTH_LOW, PEPTIDE_LENGTH_HIGH + 1):
    runners.run_template_peptides(peptide_length=i)
