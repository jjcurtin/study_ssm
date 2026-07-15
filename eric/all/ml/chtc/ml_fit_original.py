import numpy as np
import pandas as pd
import sys,time
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import StratifiedKFold,GridSearchCV, RandomizedSearchCV
from sklearn.utils import shuffle
from xgboost import XGBClassifier
from scipy.stats import randint

process_id = int(sys.argv[1])
cluster_id = int(sys.argv[2])
subid = int(sys.argv[3])
raw_feats_labels = pd.read_csv("ml_labels_alt.csv")

use_XGB = True
print(use_XGB)
windows = [0,1,3,7]
print(raw_feats_labels.shape)
lr_output = raw_feats_labels.query("subid==@subid").copy(deep=True)[['subid','day']]
xgb_output = raw_feats_labels.query("subid==@subid").copy(deep=True)[['subid','day']]
for window in windows:
    pred_str = 'w'+str(window)+'_pred'
    act_str = 'w'+str(window)+'_act'
    lr_output[pred_str] = np.nan
    xgb_output[pred_str] = np.nan
    lr_output[act_str] = np.nan
    xgb_output[act_str] = np.nan

# Hyperparameter search grids
lr_hyperparameter_grid = {'penalty':['l1','l2'],'C': [0.001,0.01,0.1,1,10,100,1000]}
xgb_hyperparameter_grid = {'n_estimators':[300,400,500,600,700], 'colsample_bytree':[0.5,1],
          'max_depth':[2,3,4,5],'min_child_weight':[2,3,4],'eta':[0.1]}
xgb_hyperparameter_grid_alt = {'n_estimators':np.arange(40,130,10),'eta':[0.1], 'colsample_bytree':[.5,1],
          'max_depth':[1,2,3],'min_child_weight':[1,2,3]}
gb_hyperparameter_grid = {'n_estimators':randint(10,500),'max_features':['sqrt', None],
          'max_depth':randint(1,10),'min_samples_leaf':randint(1,10)}

n_folds = 10
# Loop over the windows and fit models

for window in windows:
    print(window)
    label_str = "lapse_w"+str(window)
    temp = raw_feats_labels.query("subid==@subid")
    idxs = temp[temp[label_str].notnull()].index.to_numpy()
    train_feats_labels = shuffle(raw_feats_labels.query("subid != @subid").copy(deep=True).dropna(subset=label_str),random_state=110)
    #print(train_feats_labels.shape)
    test_feats_labels = raw_feats_labels.query("subid == @subid").copy(deep=True).dropna(subset=label_str)
    #print(test_feats_labels.shape)
    # Training
    train_feats = train_feats_labels.drop(labels=['day','subid','lapse_w0','lapse_w1','lapse_w3','lapse_w7'],axis=1)
    #print(train_feats.shape)
    test_feats = test_feats_labels.drop(labels=['day','subid','lapse_w0','lapse_w1','lapse_w3','lapse_w7'],axis=1)
    #print(test_feats.shape)
    # Testing 
    train_labels = train_feats_labels[label_str]
    #print(train_labels.shape)
    test_labels = test_feats_labels[label_str]
    #print(test_labels.shape)
    LRmodel = LogisticRegression(solver='liblinear')
    clfLR = GridSearchCV(LRmodel,
                param_grid=lr_hyperparameter_grid,
                cv=n_folds,n_jobs=6, scoring ='f1')
    #clfLR.fit(train_feats,train_labels)
    clfLR.fit(train_feats.drop(labels='day_of_week_num_0',axis=1),train_labels)

    # XGBmodel = XGBClassifier(n_jobs=None)
    # clfXGB = RandomizedSearchCV(XGBmodel,
    #     param_distributions=xgb_hyperparameter_grid,
    #     cv=folds, n_iter=500,
    #     scoring ='roc_auc',n_jobs=1)
    # clfXGB.fit(train_feats,train_labels)
    start_time = time.time()
    if(use_XGB):
        XGBmodel = XGBClassifier(n_jobs=1,objective='binary:logistic',validate_parameters=True)
        # clfXGB = RandomizedSearchCV(XGBmodel,
        #     param_distributions=xgb_hyperparameter_grid,
        #     cv=folds, n_iter=750, n_jobs=6,
        #     scoring ='f1',verbose=2)
        if(label_str=='lapse_w0'):
            current_grid = xgb_hyperparameter_grid_alt
        else:
            current_grid = xgb_hyperparameter_grid_alt
        clfXGB = GridSearchCV(XGBmodel,
                              param_grid=current_grid,
                              cv = n_folds,n_jobs=6,scoring='f1',verbose=1)
        clfXGB.fit(train_feats,train_labels)
        print(clfXGB.best_params_)
    else:
        GBmodel = GradientBoostingClassifier()
        clfXGB = RandomizedSearchCV(GBmodel,
            param_distributions=gb_hyperparameter_grid,
            cv=n_folds, n_iter=100, n_jobs = 4,
            scoring ='f1',verbose=2)
        clfXGB.fit(train_feats,train_labels)
    print(time.time()-start_time)
    pred_str = 'w'+str(window)+'_pred'
    act_str = 'w'+str(window)+'_act'
    #print(test_feats)
    if test_labels.shape[0]==0:
        continue
    #lr_output.loc[idxs,pred_str] = clfLR.best_estimator_.predict_proba(test_feats)[:,1]
    lr_output.loc[idxs,pred_str] = clfLR.best_estimator_.predict_proba(test_feats.drop(labels='day_of_week_num_0',axis=1))[:,1]
    lr_output.loc[idxs,act_str] = test_labels
    xgb_output.loc[idxs,pred_str] = clfXGB.best_estimator_.predict_proba(test_feats)[:,1]
    xgb_output.loc[idxs,act_str] = test_labels

lr_output['fit_type'] = 'LR'
if use_XGB:
    xgb_output['fit_type'] ='XGB'
else:
    xgb_output['fit_type'] = 'GB'

out = pd.concat([lr_output,xgb_output],ignore_index=True)
out_str = "outputs/{}_{}.csv".format(str(cluster_id),str(process_id))
out.to_csv(out_str,index=False)