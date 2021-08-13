"""
This is a script for basic processing in Perseus pipeline
"""
import numpy as np


# I Helper Functions
def get_dict_groups(df=None, lfq_str="log2 LFQ", groups=None):
    """Get dict with groups from df based on lfq_str and given groups"""
    dict_col_group = {}
    for col in list(df):
        if lfq_str in col:
            for group in groups:
                if group in col:
                    dict_col_group[col] = group
    return dict_col_group


def pre_filter(df=None, list_filter_col=None):
    """Remove row if nan for filtering columns"""
    if list_filter_col is None:
        list_filter_col = ["Only identified by site", "Reverse", "Contaminant"]
    for filter_col in list_filter_col:
        if filter_col in list(df):
            df = df[df[filter_col].isna()]
    return df


# II Main Functions
class PerseusBase:
    """Base class for Perseus Pipeline"""

    def __init__(self, df=None, dict_col_group=None, col_acc="Protein ID", col_genes="Gene Names",
                 pre_filtered=False):
        dict_group_cols = {dict_col_group[key]: [] for key in dict_col_group}
        for key in dict_col_group:
            dict_group_cols[dict_col_group[key]].append(key)
        # Dict with group associations
        self.dict_col_group = dict_col_group
        self.dict_group_cols = dict_group_cols
        self.list_groups = list(self.dict_group_cols.keys())
        # list cols
        self.list_col_lfq = list(self.dict_col_group.keys())
        list_col = [col_acc, col_genes] + self.list_col_lfq
        # Pre filter df
        if pre_filtered:
            df = pre_filter(df=df)
        df_modified = df[list_col].copy()
        df_modified.rename({col_acc: "ACC", col_genes: "Gene_Name"}, axis=1, inplace=True)
        self._df = df_modified

    def get_df_lfq(self, log2_in=True, log2_out=True):
        """Get df with just LFQ values in log2 or normal scale"""
        df_lfq = self._df[self.list_col_lfq]
        if log2_out:
            if not log2_in:
                self.list_col_lfq = ["log2 {}".format(x) for x in self.list_col_lfq]
                df_lfq = np.log2(df_lfq.replace(0, np.nan))
        else:
            if log2_in:
                self.list_col_lfq = [x.replace("log2 ", "") for x in self.list_col_lfq]
                df_lfq = 2 ** df_lfq
        df_lfq.columns = self.list_col_lfq
        return df_lfq

    def add_acc_gene(self, df=None):
        """Add UniProt accession number and gene name to df based on index"""
        df_acc_gene = self._df[["ACC", "Gene_Name"]]
        df_out = df_acc_gene.join(df, how="right")
        return df_out

