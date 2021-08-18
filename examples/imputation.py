"""
This is a script for ...
"""
import time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import examples._utils as ut
from perseuspy import PerseusPipeline, get_dict_groups

# Settings
pd.set_option('expand_frame_repr', False)  # Single line print for pd.Dataframe
FOLDER_IN = ut.FOLDER_DATA + "Liu_Dongre_2020" + ut.SEP

# I Helper Functions


# II Main Functions
def imputation1():
    """"""
    file = "PXD000485.tsv"
    df = pd.read_csv(FOLDER_IN + file)
    dict_col_group = get_dict_groups(df=df, lfq_str="LFQ intensity", groups=["OR", "PD"])
    pp = PerseusPipeline(dict_col_group=dict_col_group, df=df, col_acc="Protein IDs", col_genes="Gene names")
    df_ratio_pval = pp.run(log2_in=False).dropna()
    print(df_ratio_pval)
    pp.volcano_plot(df_ratio_pval=df_ratio_pval, col_pval="-log10 p value (OR/PD)",
                    col_ratio="log2 ratio (OR/PD)", minor_ticks=True)
    plt.show()

# III Test/Caller Functions


# IV Main
def main():
    t0 = time.time()
    imputation1()
    t1 = time.time()
    print("Time:", t1 - t0)


if __name__ == "__main__":
    main()
