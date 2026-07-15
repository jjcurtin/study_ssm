#!/bin/zsh

# Tarball the relevant files (-C removes the path in the tarball)
tar -czf input_files.tar.gz \
-C /Users/eric/repos/aud/data mle_coef_fits_lapse_not_fitted.csv day_labels.csv subject_info.csv ssm_fold_reference.json \
-C /Users/eric/repos/aud map_fit_driver.R rolling_map_fit_driver.R aud_helper_functions.R 

scp input_files.tar.gz pulick@ap2002.chtc.wisc.edu:/home/pulick
scp aud_driver.sh pulick@ap2002.chtc.wisc.edu:/home/pulick
scp aud.sub pulick@ap2002.chtc.wisc.edu:/home/pulick
scp aud_driver_rolling.sh pulick@ap2002.chtc.wisc.edu:/home/pulick
scp aud_rolling.sub pulick@ap2002.chtc.wisc.edu:/home/pulick
scp subid_list.txt pulick@ap2002.chtc.wisc.edu:/home/pulick
scp subid_list_data_eff.txt pulick@ap2002.chtc.wisc.edu:/home/pulick

#rm input_files.tar.gz