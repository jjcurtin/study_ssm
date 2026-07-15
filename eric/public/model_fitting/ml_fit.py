import numpy as np
import pandas as pd
import sys,time,json
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import StratifiedKFold,GroupKFold,GridSearchCV,StratifiedGroupKFold, RandomizedSearchCV
from sklearn.utils import shuffle
from xgboost import XGBClassifier
from scipy.stats import randint
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
run_CHTC = False

# Separate for cluster runs and local testing
# Run on cluster (CHTC)
if run_CHTC:
    process_id = int(sys.argv[1])
    cluster_id = int(sys.argv[2])
    key = sys.argv[3]
    with open('ml_fold_reference.json') as f:
        fold_reference = json.load(f)
    raw_feats_labels = pd.read_csv("ml_labels_full_width.csv")
    test_feats_labels = pd.read_csv('ml_labels_limited_width.csv')
    raw_info = pd.read_csv("subject_info.csv")
# Local testing
else:
    process_id = 1
    cluster_id = 1
    #subid = 5
    key = 'split5_rs1_fold0_id0'
    with open('../data/ml_fold_reference.json') as f:
        fold_reference = json.load(f)
    raw_feats_labels = pd.read_csv("../data/ml_labels_full_width.csv")
    test_feats_labels = pd.read_csv("../data/ml_labels_limited_width.csv")
    raw_info = pd.read_csv("../data/subject_info.csv")

current_fold_info = fold_reference[key]
train_ids = current_fold_info['train']
test_ids = current_fold_info['test']
windows = [0,1,3,7]
print("raw_feats_labels shape",raw_feats_labels.shape)

# Hyperparameter search grids
lr_hyperparameter_grid = {'penalty':['l1','l2'],'C': [0.001,0.01,0.1,1,10,100,1000]}
xgb_hyperparameter_grid = {'n_estimators':[25,50,75,100,125,150,175,200,225,250,275,300,325,350,375,400],'eta':[0.1],
                           'colsample_bytree':[.75], 'subsample':[0.5],'max_depth':[1,2,3,4,5],'min_child_weight':[1,2,3,4,5]}

n_folds = 5
job_count = 8
# Loop over the windows and fit models
results_list = []

for window in windows:
    # Process the iteration window (using the corresponding string for column identification)
    print("Window:",window)
    label_str = "lapse_w"+str(window)

    # Query down the training and testing sets based on the inputted lists
    train_set = raw_feats_labels.query("subid in @train_ids").copy(deep=True)
    test_set = test_feats_labels.query("subid in @test_ids").copy(deep=True)

    # Create the training and testing label sets, further filtered down by the predictions that are relevant for this window
    # For instance, predictions can't be made when the window extends beyond the subject's study period
    train_feats_labels = train_set.copy(deep=True).dropna(subset=label_str)
    print("train_feats_labels shape", train_feats_labels.shape)
    #print(train_feats_labels.shape)
    test_feats_labels = test_set.copy(deep=True).dropna(subset=label_str)
    print("test_feats_labels shape",test_feats_labels.shape)
    #print(test_feats_labels.shape)

    # Make the output dataframes from the filtered dataframes
    xgb_output = test_feats_labels.copy(deep=True)[['subid','day','train_width']]
    lr_output = test_feats_labels.copy(deep=True)[['subid','day','train_width']]
    xgb_output['key']=key
    lr_output['key']=key
    pred_str = 'pred'
    act_str = 'act'
    lr_output['window'] = window
    xgb_output['window'] = window

    # Training
    train_groups = train_feats_labels['subid']
    train_feats = train_feats_labels.drop(labels=['day','subid','lapse_w0','lapse_w1','lapse_w3','lapse_w7','train_width'],axis=1)
    print(train_feats.shape)
    #print(train_feats.shape)
    test_feats = test_feats_labels.drop(labels=['day','subid','lapse_w0','lapse_w1','lapse_w3','lapse_w7','train_width'],axis=1)
    #print(test_feats.shape)
    # Testing 
    train_labels = train_feats_labels[label_str]
    #print(train_labels.shape)
    test_labels = test_feats_labels[label_str]
    #print(test_labels.shape)
    # CV folds
    cv_generator = StratifiedGroupKFold(n_splits=n_folds, shuffle=True, random_state=15)
    LRmodel = LogisticRegression(solver='liblinear')
    clfLR = GridSearchCV(LRmodel,
                param_grid=lr_hyperparameter_grid,
                cv=cv_generator,n_jobs=job_count, scoring ='f1')

    clfLR.fit(train_feats.drop(labels='day_of_week_num_0',axis=1),train_labels,groups = train_groups)

    start_time = time.time()
    XGBmodel = XGBClassifier(n_jobs=1,objective='binary:logistic',validate_parameters=True)
    current_grid = xgb_hyperparameter_grid
    clfXGB = GridSearchCV(XGBmodel,
                        param_grid=current_grid,
                        cv = cv_generator,n_jobs=job_count,scoring='f1',verbose=1)
    clfXGB.fit(train_feats,train_labels,groups=train_groups)
    print(clfXGB.best_params_)
    print(time.time()-start_time)

    lr_output[pred_str] = clfLR.best_estimator_.predict_proba(test_feats.drop(labels='day_of_week_num_0',axis=1))[:,1]
    lr_output[act_str] = test_labels
    xgb_output[pred_str] = clfXGB.best_estimator_.predict_proba(test_feats)[:,1]
    xgb_output[act_str] = test_labels
    
    lr_output['fit_type'] = 'LR'
    xgb_output['fit_type'] ='XGB'

    results_list.append(lr_output)
    results_list.append(xgb_output)

out = pd.concat(results_list,ignore_index=True)
if run_CHTC:
    out_str = "outputs/{}_{}_{}.csv".format(str(cluster_id),str(process_id),str(key))
else:
    out_str = "../outputs/{}_{}.csv".format(str(cluster_id),str(process_id))
out.to_csv(out_str,index=False)