"""
This is a script for statistical tests in Perseus pipeline
"""
import time
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
import functools
import warnings


from perseuspy.per_base import PerseusBase


# I Helper Functions
def ignore_warning(simplefilter=True, category=RuntimeWarning):
    """Ignore user warning just for function"""
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            # Do something before
            if simplefilter:
                warnings.simplefilter("ignore", UserWarning)
            else:
                warnings.filterwarnings("ignore", category=category)
            # Func
            value = func(*args, **kwargs)
            # Do something after
            if simplefilter:
                warnings.simplefilter("default", UserWarning)
            else:
                warnings.filterwarnings("default", category=category)
            return value
        return wrapper
    return decorator


def _check_p_correction(method=None):
    """Check p value correction methods"""
    p_corrections = ["bonferroni", "sidak", "holm", "hommel", "fdr_bh"]
    if method is not None and method not in p_corrections:
        raise ValueError("P value correction should be one of following: " + str(p_corrections))


def _correct_p_val(p_vals=None, method=None):
    """Correct p values with given methods"""
    _check_p_correction(method=method)
    if method is not None:
        p_vals = [p if str(p) != "nan" else 1 for p in p_vals]
        cor_p_vals = multipletests(p_vals, method=method)[1]
        p_vals = [p if str(p_vals[i]) != "nan" else np.nan for i, p in enumerate(cor_p_vals)]
    return p_vals


@ignore_warning(simplefilter=False, category=RuntimeWarning)
def _ttest_no_warning(df_group_a, df_group_b):
    """Call ttest ignoring RuneTime Warning"""
    p_vals = ttest_ind(df_group_a, df_group_b, axis=1, nan_policy="omit")[1]
    return p_vals


# II Main Functions
class PerseusTests(PerseusBase):
    """Class for Perseus analysis"""
    def __init__(self, **kwargs):
        PerseusBase.__init__(self, **kwargs)

    def ttest(self, df_lfq=None, method=None, nan_policy="omit", log10_out=True):
        """Pairwise t test for groups of data frame
        In: a) df_lfq: df with lfq values (in log2 scale with values for each sample)
            b) method: Correction method for ttest {None, "bonferroni", "sidak", "holm", "hommel", "fdr_bh"}
            c) nan_policy: NaN handling {'propagate', 'raise', 'omit'}
                'propagate': returns nan if any value NaN
                'raise': throws an error
                'omit': performs the calculations ignoring nan values
            d) log10_out: Boolean to decide whether p value should be in -log10 or normal scale
        Out:a) df_pval: df with p value for each group comparison
        """
        _check_p_correction(method=method)
        pval_str = "p value "
        dict_p_vals = {}
        for a in self.list_groups:
            for b in self.list_groups:
                df_group_a = df_lfq[self.dict_group_cols[a]]
                df_group_b = df_lfq[self.dict_group_cols[b]]
                if nan_policy == "omit":
                    p_vals = _ttest_no_warning(df_group_a, df_group_b)
                else:
                    p_vals = ttest_ind(df_group_a, df_group_b, axis=1, nan_policy=nan_policy)[1]
                p_vals = _correct_p_val(p_vals=p_vals, method=method)
                dict_p_vals[pval_str + "({}/{})".format(b, a)] = p_vals
        df_pval = pd.DataFrame(dict_p_vals)
        if log10_out:
            cols = ["-log10 {}".format(x) for x in list(df_pval)]
            df_pval = -np.log10(df_pval)
            df_pval.columns = cols
        return df_pval

    def fdr_permutation(self):
        """Permutation based FDR correction

        References
        ----------
        [1] Tusher, G. V., Tibshirani, R. & Gilbert C. Significance analysis of microarrays applied to
        the ionizing radiation response. PNAS (2001)

        [2] Tyanova, S. & Cox, J. Perseus: A Bioinformatics Platform for Integrative Analysis of Proteomics
        Data in Cancer Research. Springer Protocols - Cancer Systems Biology (2010)
        """
        # TODO implement


    def anova(self):
        """"""
        # TODO implement
        pass
