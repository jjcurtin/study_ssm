import numpy as np
import pandas as pd
import sys,time
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import StratifiedKFold,GroupKFold,GridSearchCV,StratifiedGroupKFold, RandomizedSearchCV
from sklearn.utils import shuffle
from xgboost import XGBClassifier
from scipy.stats import randint
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
run_CHTC = True

if run_CHTC:
    process_id = int(sys.argv[1])
    cluster_id = int(sys.argv[2])
    subid = int(sys.argv[3])
    raw_feats_labels = pd.read_csv("ml_labels.csv")
    raw_info = pd.read_csv("subject_info.csv")
else:
    process_id = 1
    cluster_id = 1
    subid = 5
    raw_feats_labels = pd.read_csv("/Users/eric/repos/aud/data/ml_labels.csv")
    raw_info = pd.read_csv("/Users/eric/repos/aud/data/subject_info.csv")

subject_info = raw_info.query("subid==@subid")
first_ema_day = subject_info.first_morning_ema_day.iloc[0]
last_ema_day = subject_info.last_morning_ema_day.iloc[0]

use_XGB = True
print(use_XGB)
windows = [0,1,3,7]
print(raw_feats_labels.shape)

# Hyperparameter search grids
lr_hyperparameter_grid = {'penalty':['l1','l2'],'C': [0.001,0.01,0.1,1,10,100,1000]}
# xgb_hyperparameter_grid = {'n_estimators':[50,100,150,200,250],'eta':[0.2], 'colsample_bytree':[.75], 'subsample':[0.5],
#           'max_depth':[1,2,3],'min_child_weight':[1,3]}
xgb_hyperparameter_grid = {'n_estimators':[25,50,75,100,125,150,175,200,225,250, 275,300],'eta':[0.1], 'colsample_bytree':[.75], 'subsample':[0.5],
          'max_depth':[1,2,3,4],'min_child_weight':[1,3,5]}
gb_hyperparameter_grid = {'n_estimators':randint(10,500),'max_features':['sqrt', None],
          'max_depth':randint(1,10),'min_samples_leaf':randint(1,10)}

n_folds = 20
job_count = 8
# Loop over the windows and fit models
results_list = []


for window in windows:
    print(window)
    label_str = "lapse_w"+str(window)
    temp = raw_feats_labels.query("subid==@subid")
    idxs = temp[temp[label_str].notnull()].index.to_numpy()
    non_nan_days = temp.loc[idxs].day.unique()

    start_day = max(min(non_nan_days),first_ema_day)
    stop_day = min(max(non_nan_days),last_ema_day)
    non_nan_days = non_nan_days[non_nan_days >= start_day]
    non_nan_days = non_nan_days[non_nan_days <= stop_day]
    for iter_day in non_nan_days[0:1]:
        iter_idxs = raw_feats_labels.query("subid==@subid & day in @non_nan_days & day>=@iter_day").index.to_numpy()
        xgb_output = raw_feats_labels.loc[iter_idxs].copy(deep=True)[['subid','day']]
        lr_output = raw_feats_labels.loc[iter_idxs].copy(deep=True)[['subid','day']]
        pred_str = 'pred'
        act_str = 'act'
        lr_output[pred_str] = np.nan
        xgb_output[pred_str] = np.nan
        lr_output[act_str] = np.nan
        xgb_output[act_str] = np.nan
        lr_output["train_horizon"] = iter_day-1
        xgb_output["train_horizon"] = iter_day-1
        lr_output['window'] = window
        xgb_output['window'] = window

        #train_feats_labels = shuffle(raw_feats_labels.drop(iter_idxs,axis=0).copy(deep=True).dropna(subset=label_str),random_state=15)
        train_feats_labels = raw_feats_labels.drop(iter_idxs,axis=0).copy(deep=True).dropna(subset=label_str)
        print(train_feats_labels.shape)
        #print(train_feats_labels.shape)
        test_feats_labels = raw_feats_labels.loc[iter_idxs].copy(deep=True).dropna(subset=label_str)
        print(test_feats_labels.shape)
        #print(test_feats_labels.shape)
        # Training
        train_groups = train_feats_labels['subid']
        train_feats = train_feats_labels.drop(labels=['day','subid','lapse_w0','lapse_w1','lapse_w3','lapse_w7'],axis=1)
        print(train_feats.shape)
        #print(train_feats.shape)
        test_feats = test_feats_labels.drop(labels=['day','subid','lapse_w0','lapse_w1','lapse_w3','lapse_w7'],axis=1)
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
        if(use_XGB):
            XGBmodel = XGBClassifier(n_jobs=1,objective='binary:logistic',validate_parameters=True)
            current_grid = xgb_hyperparameter_grid
            clfXGB = GridSearchCV(XGBmodel,
                                param_grid=current_grid,
                                cv = cv_generator,n_jobs=job_count,scoring='f1',verbose=1)
            clfXGB.fit(train_feats,train_labels,groups=train_groups)
            print(clfXGB.best_params_)
        else:
            GBmodel = GradientBoostingClassifier()
            clfXGB = RandomizedSearchCV(GBmodel,
                param_distributions=gb_hyperparameter_grid,
                cv=n_folds, n_iter=100, n_jobs = 4,
                scoring ='f1',verbose=2)
            clfXGB.fit(train_feats,train_labels)
        print(time.time()-start_time)
        pred_str = 'pred'
        act_str = 'act'

        if test_labels.shape[0]==0:
            continue
        lr_output.loc[iter_idxs,pred_str] = clfLR.best_estimator_.predict_proba(test_feats.drop(labels='day_of_week_num_0',axis=1))[:,1]
        lr_output.loc[iter_idxs,act_str] = test_labels
        xgb_output.loc[iter_idxs,pred_str] = clfXGB.best_estimator_.predict_proba(test_feats)[:,1]
        xgb_output.loc[iter_idxs,act_str] = test_labels
        
        lr_output['fit_type'] = 'LR'
        if use_XGB:
            xgb_output['fit_type'] ='XGB'
        else:
            xgb_output['fit_type'] = 'GB'
        results_list.append(lr_output)
        results_list.append(xgb_output)

out = pd.concat(results_list,ignore_index=True)
if run_CHTC:
    out_str = "outputs/{}_{}.csv".format(str(cluster_id),str(process_id))
else:
    out_str = "{}_{}.csv".format(str(cluster_id),str(process_id))
out.to_csv(out_str,index=False)