

#Making an isolated virtualenv for trgt
conda create -n mtrgtpy python=3
conda activate mtrgtpy
#!/bin/bash

# Create and activate a new conda environment
conda create -n mtrgtpy python=3 -y
conda activate mtrgtpy

# Navigate to the Truvari directory and install it
cd truvari  ## available in pypr

pip install .

pip install truvari

# Navigate back to the original directory
cd -

# Install trgt and its dependencies
pip install .

# List the number of installed packages
conda list | wc -l


# manually installing truvari v4.0-dev
cd ../truvari/
python3 -m pip install .
cd -
# installing trgt and its dependencies
python3 -m pip install .
python3 -m pip freeze | wc -l