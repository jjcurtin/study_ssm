#!/bin/bash

# untar input files
tar -xzf input_ml_files_data_eff.tar.gz

# make outputs folder
mkdir outputs

# run your script
python3 ml_fit_data_eff.py $1 $2 $3