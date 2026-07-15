The auROC and auPRC values from 15 repeats of 5-fold stratified cross valdiation can be found in the files auroc_vals.csv and auprc_vals.csv, respectively.

The columns refer to the following:
window - The prediction window length (0 is same-day, 3 is within 3 days, 7 is within 7 days).
width - The width of the data availability window (e.g., 15, 30, 45, etc.)
fold_id - The id for the fold-within-repeat (i.e., a value from 0-4)
rs_id - The id for the repeat (i.e., a value from 0-14)
ssm - The auROC/auPRC associated with the state space model
lr - The auROC/auPRC associated with the logistic regression model
xgb - The auROC/auPRC associated with the XGB model