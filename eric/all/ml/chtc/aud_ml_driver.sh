#!/bin/bash

# untar input files
tar -xzf input_ml_files.tar.gz

# make outputs folder
mkdir outputs

# run your script
python3 ml_fit.py $1 $2 $3