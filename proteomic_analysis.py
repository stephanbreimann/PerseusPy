"""
This is a script for ...
"""
import os
import time
import pandas as pd
import matplotlib.pyplot as plt

from perseuspy import PerseusPipeline, get_dict_groups

# Settings
pd.set_option('expand_frame_repr', False)  # Single line print for pd.Dataframe
FOLDER_DATA = os.path.dirname(os.path.abspath(__file__)) + "/data/"

# I Helper Functions


# II Main Functions
def perseus_analysis():
    """"""
    # Load data
    df = pd.read_csv(FOLDER_DATA + "npc_symptomatic_mice.csv")
    # 1. Perseus Analysis
    groups = ["Npc1-/-_8weeks", "WT_8weeks"]
    dict_col_group = get_dict_groups(df=df, lfq_str="Log2 LFQ", groups=groups)
    pp = PerseusPipeline(df=df, dict_col_group=dict_col_group)
    df_ratio_pval = pp.run(gmean=False)
    print(df_ratio_pval)
    pp.volcano_plot(df_ratio_pval=df_ratio_pval,
                    col_pval="-log10 p value (WT_8weeks/Npc1-/-_8weeks)",
                    col_ratio="log2 ratio (WT_8weeks/Npc1-/-_8weeks)",
                    th_filter=(0.05, 0.5),
                    th_text=(0.05, -2.5, 2.5),
                    force=(0.9, 0.50, 0.25),
                    avoid_conflict=0.2,
                    precision=0.01,
                    box=True,
                    verbose=True,
                    label_bold=False,
                    label_size=10,
                    filled_circle=True,
                    title="Test_Volcano_mean",
                    fig_format="png")
    plt.show()

# III Test/Caller Functions


# IV Main
def main():
    t0 = time.time()
    perseus_analysis()
    t1 = time.time()
    print("Time:", t1 - t0)


if __name__ == "__main__":
    main()
