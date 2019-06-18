#!/bin/bash

# write name of environment to environment.yml file
echo "name: ${1}" > test.yml
tail -n +2 environment.yml >> test.yml
mv test.yml environment.yml

# activate conda shell so conda command can be used
source $(conda info --base)/etc/profile.d/conda.sh

# create new conda environment from environment.yml file
conda env create -f environment.yml

# if successful load environemnt variable setup and tear down scripts
if [ $? -eq 0 ]; then

    PROJECT_DIR="`pwd`"
    PYTHONPATH=$PYTHONPATH

    conda activate $1

    cd $CONDA_PREFIX
    mkdir -p ./etc/conda/activate.d
    mkdir -p ./etc/conda/deactivate.d

    # set up script
    cat > ./etc/conda/activate.d/env_vars.sh << EOF
    #!/bin/bash

    export PROJECT_DIR="${PROJECT_DIR}"
    export PYTHONPATH="\$PROJECT_DIR:\$PYTHONPATH"
EOF

    # tear down script
    cat > ./etc/conda/deactivate.d/env_vars.sh << EOF
    #!/bin/bash

    unset PROJECT_DIR
    export PYTHONPATH="$PYTHONPATH"
EOF

    conda deactivate
fi
