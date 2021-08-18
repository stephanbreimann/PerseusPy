"""
This is the main script of the Perseus pipeline
Ref.: Tyanova et al., 2016 (The Perseus computational platform for comprehensive analysis of (prote)omics data)
"""
import time
import pandas as pd
import numpy as np

from perseuspy.per_comput import PerseusComputations
from perseuspy.per_plots import PerseusPlots
from perseuspy.per_test import PerseusTests
import perseuspy._utils as ut


# I Helper Functions
def check_log2_scale_of_lfq(df_lfq=None, th_max_log2=100):
    """"""
    max_intensity = np.round(max(df_lfq.values.flatten()), 2)
    if max_intensity > th_max_log2:
        error = f"Maximum intensity in df ({max_intensity}) is exceeding 'th_max_log2' ({th_max_log2})." \
                f"\nValues are probably not in log2 scale. If yes, increase 'th_max_log2' to continue."
        raise ValueError(error)


# TODO heavy check input df
# II Main Functions
class PerseusPipeline(PerseusComputations, PerseusTests, PerseusPlots):    # PerseusNormalization, PerseusImputation,
    """Class for Perseus analysis"""
    def __init__(self, df=None, dict_col_group=None, col_acc=ut.COL_ACC, col_genes=ut.COL_GENE,
                 pre_filtered=False, groups=None):
        """
        Class for proteomic analysis pipeline as performed in Perseus.

        Parameters
        ----------
        df: pd.DataFrame with unique protein ids (ACC), protein/gene name and expression intensities
            (e.g., by label free quantification (lfq)) for n>=2 groups to compare with
        dict_col_group: dict with intensity column to group names
        col_acc: {str} default "Protein ID".
            Column name in df for unique protein identifier
        col_genes: {str} default "Gene Names".
            Column name in df for gene names
        pre_filtered: {bool} default False. Specify whether values in df are already filtered for
            "Only identified by site", "Reverse", "Contaminant". If False, df will be filtered.
        groups: {list} list with group names {strings}
        """
        kwargs = dict(df=df,
                      dict_col_group=dict_col_group,
                      col_acc=col_acc,
                      col_genes=col_genes,
                      pre_filtered=pre_filtered,
                      groups=groups)
        PerseusComputations.__init__(self, **kwargs)
        PerseusTests.__init__(self, **kwargs)
        PerseusPlots.__init__(self, **kwargs)

    def run(self, log2_in=True, log2_max=100):
        """Run perseuspy pipeline to get df_ratio_pval:
            df_lfq -> df_lfq_mean -> df_ratio + df_pval

        Parameters
        ----------
        log2_in: {bool} True. Specify whether intensity values in df are log2 transformed or not.
        log2_max: {int} default 100. Maximum value to decide if values are log scaled or normal scaled
        """
        # 1.1 LFQ Processing (df_lfq -> df_ratio)

        df_lfq = self.get_df_lfq(log2_in=log2_in)
        check_log2_scale_of_lfq(df_lfq=df_lfq, th_max_log2=log2_max)
        df_lfq_mean = self.get_df_lfq_mean(df_lfq=df_lfq, remove_nan=False)
        df_ratio = self.get_df_ratio(df_lfq_mean=df_lfq_mean)
        # 1.2 Statistical tests (df_lfq -> df_pval)
        df_pval = self.ttest(df_lfq=df_lfq)
        # 1.3 Join ratio and statistical analysis
        df_ratio_pval = df_ratio.join(df_pval)
        df_ratio_pval = self.add_acc_gene(df_ratio_pval)
        return df_ratio_pval

