#!/bin/zsh

# Tarball the relevant files (-C removes the path in the tarball)
tar -czf input_ml_files_data_eff.tar.gz \
-C /Users/eric/repos/aud/data ml_labels_full_width.csv ml_labels_limited_width.csv subject_info.csv \
-C /Users/eric/repos/aud/ml/chtc ml_fit_data_eff.py ml_fold_reference.json

scp input_ml_files_data_eff.tar.gz pulick@ap2002.chtc.wisc.edu:/home/pulick
scp aud_ml_driver_data_eff.sh pulick@ap2002.chtc.wisc.edu:/home/pulick
scp aud_ml_data_eff.sub pulick@ap2002.chtc.wisc.edu:/home/pulick
scp subid_ml_list_data_eff.txt pulick@ap2002.chtc.wisc.edu:/home/pulick

#rm input_files.tar.gz