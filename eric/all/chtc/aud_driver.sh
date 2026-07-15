#!/bin/bash

# untar input files
tar -xzf input_files.tar.gz

# make outputs folder
mkdir outputs

# run your script
Rscript map_fit_driver.R $1 $2 $3 $4