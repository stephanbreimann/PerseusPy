"""
This is a script for basic processing in Perseus pipeline
"""
import numpy as np
import warnings

import perseuspy._utils as ut


# I Helper Functions
def _check_lfq_str(df=None, lfq_str=None, groups=None):
    """"""
    cols_lfq = []
    for col in list(df):
        if lfq_str in col:
            for group in groups:
                if group in col:
                    cols_lfq.append(col)
    if len(cols_lfq) < len(groups):
        raise ValueError("'lfq_str'({}) does not match with given data".format(lfq_str))
    if len(cols_lfq) == len(groups):
        warnings.warn("'lfq_str'({}) might not match with sample columns: \n"
                      "{}".format(lfq_str, cols_lfq))


def _check_df(df=None, col_acc=None, col_genes=None):
    """"""
    if col_acc not in list(df):
        raise ValueError("'col_acc' ({}) not in given data: {}".format(col_acc, list(df)))
    if col_genes not in list(df):
        raise ValueError("'col_genes' ({}) not in given data: {}".format(col_genes, list(df)))


def _pre_filter(df=None, list_filter_col=None):
    """Remove row if nan for filtering columns"""
    if list_filter_col is None:
        list_filter_col = ["Only identified by site", "Reverse", "Contaminant"]
    for filter_col in list_filter_col:
        if filter_col in list(df):
            df = df[df[filter_col].isna()]
    return df


# II Main Functions
def get_dict_groups(df=None, lfq_str=ut.STR_LOG2_INTENSITY, groups=None):
    """Get dict with groups from df based on lfq_str and given groups"""
    _check_lfq_str(df=df, lfq_str=lfq_str, groups=groups)
    dict_col_group = {}
    for col in list(df):
        if lfq_str in col:
            for group in groups:
                if group in col:
                    dict_col_group[col] = group
    return dict_col_group


class PerseusBase:
    """Base class for Perseus Pipeline"""

    def __init__(self, df=None, dict_col_group=None, col_acc="Protein ID", col_genes="Gene Names",
                 pre_filtered=False, groups=None):
        _check_df(df=df, col_acc=col_acc, col_genes=col_genes)
        dict_group_cols = {dict_col_group[key]: [] for key in dict_col_group}
        for key in dict_col_group:
            dict_group_cols[dict_col_group[key]].append(key)
        # Dict with group associations
        self.dict_col_group = dict_col_group
        self.dict_group_cols = dict_group_cols
        if groups is None:
            self.list_groups = list(self.dict_group_cols.keys())
        else:
            self.list_groups = groups
        # list cols
        self.list_col_lfq = list(self.dict_col_group.keys())
        list_col = [col_acc, col_genes] + self.list_col_lfq
        # Pre filter df
        if pre_filtered:
            df = _pre_filter(df=df)
        df_modified = df[list_col].copy()
        df_modified.rename({col_acc: "ACC", col_genes: "Gene_Name"}, axis=1, inplace=True)
        self._df = df_modified

    def get_df_lfq(self, log2_in=True, log2_out=True):
        """Get df with just LFQ values in log2 or normal scale"""
        df_lfq = self._df[self.list_col_lfq]
        adjust = False
        if log2_out:
            if not log2_in:
                f = lambda x: "log2 {}".format(x)
                df_lfq = np.log2(df_lfq.replace(0, np.nan))
                adjust = True
        else:
            if log2_in:
                f = lambda x: x.replace("log2 ", "")
                df_lfq = 2 ** df_lfq
                adjust = True
        if adjust:
            self.list_col_lfq = [f(x) for x in self.list_col_lfq]
            self.dict_col_group = {col: f(self.dict_col_group[col]) for col in self.dict_col_group}
            self.dict_group_cols = {group: [f(x) for x in self.dict_group_cols[group]]
                                    for group in self.dict_group_cols}
        df_lfq.columns = self.list_col_lfq
        return df_lfq

    def add_acc_gene(self, df=None):
        """Add UniProt accession number and gene name to df based on index"""
        df_acc_gene = self._df[["ACC", "Gene_Name"]]
        df_out = df_acc_gene.join(df, how="right")
        return df_out

