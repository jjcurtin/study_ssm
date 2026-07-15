#!/bin/zsh

# Tarball the relevant files (-C removes the path in the tarball)
tar -czf input_ml_files.tar.gz \
-C /Users/eric/repos/aud/data ml_labels.csv subject_info.csv \
-C /Users/eric/repos/aud/ml/chtc ml_fit.py

scp input_ml_files.tar.gz pulick@ap2002.chtc.wisc.edu:/home/pulick
scp aud_ml_driver.sh pulick@ap2002.chtc.wisc.edu:/home/pulick
scp aud_ml.sub pulick@ap2002.chtc.wisc.edu:/home/pulick
scp subid_ml_list.txt pulick@ap2002.chtc.wisc.edu:/home/pulick

#rm input_files.tar.gz