"""
This is a script for computations of Perseus pipeline
"""
import math
import numpy as np
import pandas as pd
from scipy import stats

import perseuspy._utils as ut
from perseuspy.per_base import PerseusBase


# I Helper Functions
def _inverse_log(values=None, base=2):
    """Inverse logarithmize values for given base"""
    # 2^y=x, log2(x) = y, where x is ratio and y the log2 ratio
    de_log_val = [base ** val for val in values]
    return de_log_val


def _log(values=None, base=2):
    """Calculate log values"""
    log_val = [math.log(val, base) for val in values]
    return log_val


def _ratio(group_a, group_b, log2_in=True, log2_out=True):
    """Calculate ratio for log2 transformed data"""
    if log2_in:
        ratio_log2 = group_a - group_b  # log(a/b) = log(a) - log(b)
        if log2_out:
            return ratio_log2
        else:
            return _inverse_log(values=ratio_log2, base=2)
    else:
        ratio = group_a / group_b
        if log2_out:

            return _log(values=ratio, base=2)
        else:
            return ratio


# II Main Functions
class PerseusComputations(PerseusBase):
    """Class for Perseus analysis"""
    def __init__(self, **kwargs):
        PerseusBase.__init__(self, **kwargs)

    def get_df_lfq_mean(self, df_lfq=None, log2_in=True, remove_nan=False):
        """Calculation of mean of groups for label free quantification (LFQ)
        In: a) df_lfq: df with lfq values
            b) log2_lfq: boolean to indicate if lfq values are in log2 scale
            d) remove_nan: boolean to indicate if nan should be removed if any group is completely missing
        Out:a) df_lfq_mean: df with mean lfq values for each group"""
        if log2_in:
            lfq_str = ut.STR_LOG2_INTENSITY
        else:
            lfq_str = ut.STR_INTENSITY
        dict_avg = {}
        for group in self.dict_group_cols:
            group_col = self.dict_group_cols[group]
            print(self.dict_group_cols)
            df_group = df_lfq[group_col]
            # TODO invert log2 values
            # Calculate mean without 0
            group_mean = df_group.replace(0, np.nan).mean(axis=1, skipna=True).replace(np.nan, 0)
            dict_avg["{} {}".format(lfq_str, group)] = group_mean
        df_lfq_mean = pd.DataFrame(dict_avg)
        df_lfq_mean = df_lfq_mean.replace(0, np.nan)
        if remove_nan:
            df_lfq_mean = df_lfq_mean[~df_lfq_mean.isna().any(axis=1)]
        return df_lfq_mean

    def get_df_ratio(self, df_lfq_mean=None, log2_in=True, log2_out=True):
        """Get df with ratios for group comparison of mean lfq values
        In: a) df_lfq_mean: df with mean lfq values for each group
            b) log2_in: boolean to indicate if df_lfq_mean in log2 scale
            c) log2_out: boolean to indicate if return df in log2 scale
        Out:a) df_ratio: df with ratios for individual group comparison"""
        if log2_out:
            ratio_str = ut.STR_LOG2_RATIO
        else:
            ratio_str = ut.STR_RATIO
        dict_group_col = dict(zip(self.list_groups, list(df_lfq_mean)))
        dict_ratio = {}
        ratio_pairs = []
        for a in self.list_groups:
            for b in self.list_groups:
                if a != b and {a, b} not in ratio_pairs:
                    df_group_a = df_lfq_mean[dict_group_col[a]]
                    df_group_b = df_lfq_mean[dict_group_col[b]]
                    ratio = _ratio(df_group_a, df_group_b,
                                   log2_in=log2_in,
                                   log2_out=log2_out)
                    ratio_pairs.append({a, b})
                    dict_ratio[ratio_str + " ({}/{})".format(a, b)] = ratio
        df_ratio = pd.DataFrame(dict_ratio, index=df_lfq_mean.index)     # Set index of data frame
        return df_ratio
