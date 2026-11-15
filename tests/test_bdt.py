import numpy as np
import pytest
from higgs_review.bdt_validation import validate_arrays, evaluate_holdout


def test_reject_mass_feature():
    with pytest.raises(ValueError):
        validate_arrays(np.ones((120,1)), np.array([0,1]*60), ['pt1_norm'])


def test_holdout_reports_no_unfounded_significance():
    rng=np.random.default_rng(14)
    labels=np.array([0]*200+[1]*200)
    x=rng.normal(size=(400,2))
    x[labels==1]+=0.8
    r=evaluate_holdout(x,labels,['iso1','adphi'],seed=44)
    assert 0.5<r['holdout_auc']<=1
    assert r['test_n']==80
    assert 0<=r['holdout_background_acceptance']<=1
    assert 'significance' in r['warning']
