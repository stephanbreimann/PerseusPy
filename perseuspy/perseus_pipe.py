"""
This is the main script of the Perseus pipeline
Ref.: Tyanova et al., 2016 (The Perseus computational platform for comprehensive analysis of (prote)omics data)
"""
import time
import pandas as pd

from perseuspy.per_comput import PerseusComputations
from perseuspy.per_plots import PerseusPlots
from perseuspy.per_test import PerseusTests


# I Helper Functions


# II Main Functions
class PerseusPipeline(PerseusComputations, PerseusTests, PerseusPlots):    # PerseusNormalization, PerseusImputation,
    """Class for Perseus analysis"""
    def __init__(self, **kwargs):
        PerseusComputations.__init__(self, **kwargs)
        PerseusTests.__init__(self, **kwargs)
        PerseusPlots.__init__(self, **kwargs)

    def run(self, gmean=True):
        """Run perseuspy pipeline to get df_ratio_pval:
            df_lfq -> df_lfq_mean -> df_ratio + df_pval"""
        # 1.1 LFQ Processing (df_lfq -> df_ratio)
        df_lfq = self.get_df_lfq()
        df_lfq_mean = self.get_df_lfq_mean(df_lfq=df_lfq, remove_nan=False, gmean=gmean)
        df_ratio = self.get_df_ratio(df_lfq_mean=df_lfq_mean)
        # 1.2 Statistical tests (df_lfq -> df_pval)
        df_pval = self.ttest(df_lfq=df_lfq)
        # 1.3 Join ratio and statistical analysis
        df_ratio_pval = df_ratio.join(df_pval)
        df_ratio_pval = self.add_acc_gene(df_ratio_pval)
        return df_ratio_pval

