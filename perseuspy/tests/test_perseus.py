"""
This is a script for testing PerseusPipeline
"""
import pandas as pd
import numpy as np
import pytest

import perseuspy._utils as ut
from perseuspy import PerseusPipeline, get_dict_groups

FOLDER_IN = ut.FOLDER_DATA + "test_data" + ut.SEP


# Test input
@pytest.fixture
def df_lfq():
    return pd.read_csv(FOLDER_IN + "df_lfq.tsv")


@pytest.fixture
def df_log2_lfq():
    return pd.read_csv(FOLDER_IN + "df_log2_lfq.csv")


# Corrupted input


# Tests
def test_perseus_non_log_input(df_lfq):
    dict_col_group = get_dict_groups(df=df_lfq, lfq_str="LFQ intensity", groups=["OR", "PD"])
    pp = PerseusPipeline(dict_col_group=dict_col_group, df=df_lfq,
                         col_acc="Protein IDs", col_genes="Gene names")
    df_ratio_pval = pp.run(log2_in=True)


def test_perseus_log_input(df_log2_lfq):
    dict_col_group = get_dict_groups(df=df_log2_lfq, lfq_str="Log2 LFQ", groups=["WT", "Npc1-/-"])
    pp = PerseusPipeline(dict_col_group=dict_col_group, df=df_log2_lfq,
                         col_acc="Protein ID", col_genes="Gene Names")
    df_ratio_pval = pp.run(log2_in=False)
