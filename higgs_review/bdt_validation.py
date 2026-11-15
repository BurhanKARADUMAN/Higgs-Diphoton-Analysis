"""Holdout-validated BDT classifier diagnostics for precomputed user features.

Not a Higgs discovery statistic. Does not accept candidate mass / mass-proxy
features by default. All samples must be legitimately split and documented.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import numpy as np

DEFAULT_FEATURES = ('aeta1','aeta2','adeta','adphi','costh','iso1','iso2')
FORBIDDEN = {'mass','mgg','myy','m_gg','m_yy','m_gamma_gamma','pt1_norm','pt2_norm'}


def validate_arrays(x,y,feature_names):
    a = np.asarray(x,dtype=float)
    labels = np.asarray(y,dtype=int)
    names = tuple(feature_names)
    if a.ndim != 2 or labels.shape != (len(a),) or len(names) != a.shape[1]:
        raise ValueError('Features and labels shapes do not match')
    if len(a)<100 or np.any(~np.isfinite(a)) or set(np.unique(labels)) != {0,1} or min(np.count_nonzero(labels==0),np.count_nonzero(labels==1)) < 20:
        raise ValueError('Need >=100 finite examples and both labels 0/1')
    bad = {n for n in names if n.lower() in FORBIDDEN or 'mass' in n.lower()}
    if bad:
        raise ValueError(f'Forbidden mass / mass-proxy features in default classifier: {sorted(bad)}')
    return a, labels


def evaluate_holdout(x,y,feature_names=DEFAULT_FEATURES,*,seed=42,target_fpr=0.1):
    """Train/validation/test split, threshold on validation BKG only, test once."""
    from sklearn.ensemble import HistGradientBoostingClassifier
    from sklearn.metrics import roc_auc_score, confusion_matrix
    from sklearn.model_selection import train_test_split
    a,labels = validate_arrays(x,y,feature_names)
    if not 0 < target_fpr < 1:
        raise ValueError('target_fpr must be in (0,1)')
    xtrain,xtemp,ytrain,ytemp=train_test_split(a,labels,test_size=.4,random_state=seed,stratify=labels)
    xval,xtest,yval,ytest=train_test_split(xtemp,ytemp,test_size=.5,random_state=seed+1,stratify=ytemp)
    model=HistGradientBoostingClassifier(max_iter=120,max_leaf_nodes=15,l2_regularization=2.0,learning_rate=.07,random_state=seed)
    model.fit(xtrain,ytrain)
    val_score=model.predict_proba(xval)[:,1]
    test_score=model.predict_proba(xtest)[:,1]
    threshold=float(np.quantile(val_score[yval==0],1-target_fpr))
    pred=test_score>threshold
    tn,fp,fn,tp=confusion_matrix(ytest,pred,labels=[0,1]).ravel()
    return {
        'features':list(feature_names),'seed':seed,
        'split':'60% train / 20% validation / 20% test; stratified',
        'training_n':len(ytrain),'validation_n':len(yval),'test_n':len(ytest),
        'threshold_selected_from':'validation background quantile; NEVER test set',
        'target_validation_fpr':target_fpr,'threshold':threshold,
        'validation_auc':float(roc_auc_score(yval,val_score)),
        'holdout_auc':float(roc_auc_score(ytest,test_score)),
        'holdout_signal_efficiency':float(tp/(tp+fn)),
        'holdout_background_acceptance':float(fp/(fp+tn)),
        'warning':'Classification diagnostics only: no luminosity, MC process normalization, systematics, mass-shape decorrelation, or discovery significance.'
    }


def main(argv=None):
    p=argparse.ArgumentParser(description='BDT holdout check on your labeled feature CSV; not a Higgs significance calculation')
    p.add_argument('--csv',required=True,help='CSV with label (0 background,1 signal) and approved features')
    p.add_argument('--features',default=','.join(DEFAULT_FEATURES))
    p.add_argument('--out',default='results/bdt_holdout.json')
    p.add_argument('--seed',type=int,default=42)
    a=p.parse_args(argv)
    features=tuple(x.strip() for x in a.features.split(','))
    with open(a.csv,encoding='utf-8',newline='') as f:
        rows=list(csv.DictReader(f))
    if not rows:
        p.error('No rows in input CSV')
    try:
        x=[[float(r[n]) for n in features] for r in rows]
        y=[int(r['label']) for r in rows]
    except (KeyError,ValueError) as exc:
        p.error(f'Invalid CSV feature or label: {exc}')
    result=evaluate_holdout(x,y,features,seed=a.seed)
    output=Path(a.out);output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))

if __name__=='__main__':
    main()
