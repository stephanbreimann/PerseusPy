"""
This is a script for plotting class of Perseus pipeline
"""
import os
import time
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
import matplotlib.patches as mpatches
import seaborn as sns
import math
import plotly.express as px
from adjustText import adjust_text

import perseuspy._utils as ut


# Settings
pd.set_option('expand_frame_repr', False)  # Single line print for pd.Dataframe

COLOR_UP = "firebrick"
COLOR_DOWN = "dodgerblue"
COLOR_NOT_SIG = "gray"
COLOR_TH = "black"


# I Helper Functions
def _check_col(df, col=None):
    """Check if col in df"""
    if col not in list(df):
        raise ValueError("{} should be one of following: {}".format(col, list(df)))


def _check_gene_list(df=None, gene_list=None):
    """Check if gene in df"""
    all_genes = df["Gene_Name"].tolist()
    if gene_list is not None:
        not_in_all_genes = [x for x in gene_list if x not in all_genes]
        if len(not_in_all_genes):
            raise ValueError("Genes from 'gene_list' not in 'df_ratio_pval': {}".format(not_in_all_genes))


def _check_gene_values(gene=None, x=None, y=None):
    """Check if x or y is missing value"""
    if str(x) == "nan" or str(y) == "nan":
        raise ValueError("Missing values for gene '{}'".format(gene))


# TODO adjust for proper checking interface
def _check_log_scales(x_min=None, x_max=None):
    """"""
    if x_min < -100 or x_max > 100:
        raise ValueError("Check if values in log scale")


# Filter functions
def _color_filter(df=None, th_p=2.0, th_ratio=0.5, col_ratio=None, col_pval=None, gene_list=None):
    """Classify significant values by color"""
    colors = []
    for i, row in df.iterrows():
        gene = row["Gene_Name"]
        ratio = row[col_ratio]
        p_val = row[col_pval]
        if gene_list is not None and gene in gene_list:
            if ratio > 0:
                colors.append(COLOR_UP)
            else:
                colors.append(COLOR_DOWN)
        elif p_val >= th_p:
            if ratio >= th_ratio:
                colors.append(COLOR_UP)
            elif ratio <= -th_ratio:
                colors.append(COLOR_DOWN)
            else:
                colors.append(COLOR_NOT_SIG)
        else:
            colors.append(COLOR_NOT_SIG)
    return colors


def _label_filter(df=None, gene_list=None, th_p_text=2.0, th_neg_ratio=-0.5, th_pos_ratio=0.5,
                  col_ratio=None, col_pval=None, avoid_conflict=0.5):
    """Classify significant values by color"""
    # Adjust avoid factor
    if avoid_conflict > 1:
        avoid_conflict = avoid_conflict/100
    ac = 1 - avoid_conflict
    # Filter labels
    labels = []
    for i, row in df.iterrows():
        gene = row["Gene_Name"]
        ratio = row[col_ratio]
        p_val = row[col_pval]
        x, y = ratio, p_val
        if (gene_list is not None and gene in gene_list) or (p_val > th_p_text and (ratio < th_neg_ratio or ratio > th_pos_ratio)):
            label = (gene, x, y)
            labels.append(label)
        # Save significant but not to be shown genes to avoid overlapping conflict
        elif p_val > th_p_text*0.75 and (ratio < th_neg_ratio * ac or ratio > th_pos_ratio * ac):
            label = ("", x, y)
            labels.append(label)
    return labels


def _set_labels(labels=None, objects=None, fig_format="png", label_size=8, precision=0.01, label_bold=False,
                force_points=0.75, force_text=0.75, force_objects=0.25, box=False, alpha=0.85, verbose=True,
                th_filter=None):
    """Set labels of genes automatically
    In: a) labels: list of label tuples (label, x, y)
        b) objects: list of objects that should be considered for placement
        c) fig_format: format of plot
        d1) force_text (float): the repel force from texts is multiplied by this
            value; default (0.1, 0.25)
        d2) force_points (float): the repel force from points is multiplied by this
            value; default (0.2, 0.5)
        d3) force_objects (float): same as other forces, but for repelling
            additional objects; default (0.1, 0.25)
    """
    fontdict = dict(size=label_size)
    if label_bold:
        fontdict.update(weight="bold")
    # Helvetica, Arial
    props = dict(boxstyle='round', alpha=alpha, edgecolor="white")
    texts = []
    for label, x, y in labels:
        if abs(x) < th_filter[1] or y < th_filter[0]:
            color = COLOR_NOT_SIG
        elif x < 0:
            color = COLOR_DOWN
        else:
            color = COLOR_UP
        props.update(dict(facecolor=color))
        if label is not None:
            _check_gene_values(gene=label, x=x, y=y)
            if box:
                fontdict.update(dict(color="white"))
                texts.append(plt.text(x, y, label, fontdict=fontdict, bbox=props))
            else:
                fontdict.update(dict(color="black"))
                texts.append(plt.text(x, y, label, fontdict=fontdict))
    if verbose:
        print("{} elements have to be iteratively placed".format(len(texts)))
    adjust_text(texts,
                add_objects=objects,
                precision=precision,
                force_points=force_points,
                force_text=force_text,
                force_objects=force_objects,
                arrowprops=dict(arrowstyle="-", color='gray', lw=0.5),
                save_format=fig_format)


# II Main Functions
class PerseusPlots:
    """Class for plotting proteomics plots"""
    def __init__(self, **kwargs):
        pass

    @staticmethod
    def volcano_plot(df_ratio_pval=None, col_ratio=None, col_pval=None, gene_list=None, title=None,
                     th_filter=(0.05, 0.5), th_text=None, precision=0.01, force=(0.5, 0.5, 0.25), avoid_conflict=0.25,
                     fig_format="png", verbose=True, loc_legnd=2,
                     filled_circle=True, box=True, label_bold=False, label_size=8, minor_ticks=True):
        """Calculate p value by a two sample ttest via FDR by Benjamini Hochberg and show volcano plot
        In: a) df_ratio_pval: df with p values and ratio
            b1) col_ratio: column from df_ratio_pval to show on x-axis
            b2) col_pval: column from df_ratio_pval to show on y-axis
            c1) th_filter: tuple for filtering thresholds of p_val and ratio (p_val can be given in normal scaled)
            c1) th_text: tuple for filtering thresholds of p_val, lower_ratio, and upper_ratio
            d1) force: tuple with repulsion force (points, text, objects) to modify text
                (the higher the more distributed)
            d2) avoid_conflict: percentage of conflict to avoid
                (the higher, the more slower; 0/1 means that no resp. all additional elements are considered)
            e) title: title of plot
            f1) fig_format: format of plot if saved
            f2) show: boolean to decide whether plot should be shown or saved
        Out:a) Volcano plot saved in df_results or shown
        Notes
        -----

        """
        # TODO Adjust font Helvetica, Arial
        # TODO Include labeling for selected genes (e.g., for GO terms)
        # 1. Download all fonts to Fedora
        # dnf search microsoft windows fonts
        # https://www.reddit.com/r/Fedora/comments/e9ig9m/how_can_i_install_ms_fonts_particularly_times_new/
        # 2. Update fonts for matplotlib
        # https://scentellegher.github.io/visualization/2018/05/02/custom-fonts-matplotlib.html
        plt.rcParams["font.family"] = "sans-serif"
        plt.rcParams["font.sans-serif"] = "Arial"
        # http://physicalmodelingwithpython.blogspot.com/2015/06/making-plots-for-publication.html
        if fig_format == "pdf":
            mpl.rcParams['pdf.fonttype'] = 42
        elif "svg" in fig_format:
            mpl.rcParams['svg.fonttype'] = 'none'
        # Check and Adjust parameters
        for col in [col_ratio, col_pval]:
            _check_col(df_ratio_pval, col=col)
        _check_gene_list(df=df_ratio_pval, gene_list=gene_list)
        if title is None:
            title = "Volcano Plot for KO vs WT"
        th_p, th_ratio = th_filter
        if th_text is None:
            th_text = (th_p, -th_ratio, th_ratio)
        th_p_text, th_neg_ratio, th_pos_ratio = th_text
        # Assumed that p value is given in normal scale
        if th_p < 0.5:
            th_p = -np.log10(th_p)
        if th_p_text < 0.5:
            th_p_text = -np.log10(th_p_text)
        # Plotting settings
        kwargs_filter = dict(df=df_ratio_pval, col_ratio=col_ratio, col_pval=col_pval, gene_list=gene_list)
        colors = _color_filter(**kwargs_filter,
                               th_p=th_p,
                               th_ratio=th_ratio)
        labels = _label_filter(**kwargs_filter,
                               th_p_text=th_p_text,
                               th_neg_ratio=th_neg_ratio,
                               th_pos_ratio=th_pos_ratio,
                               avoid_conflict=avoid_conflict)
        x_min, x_max = math.floor(df_ratio_pval[col_ratio].min()) - 1, math.ceil(df_ratio_pval[col_ratio].max()) + 1
        print(x_min, x_max)
        _check_log_scales(x_min=x_min, x_max=x_max)
        y_max = 1.1 * df_ratio_pval[col_pval].max()
        # Plotting
        dict_scatter = dict(y=col_pval, x=col_ratio, kind="scatter", figsize=(5, 5))
        if filled_circle:
            dict_scatter.update(dict(color=colors))
        else:
            dict_scatter.update(dict(color="none", edgecolor=colors))
        df_ratio_pval.plot(**dict_scatter)
        if minor_ticks:
            plt.xticks(ticks=range(x_min, x_max + 1), labels=range(x_min, x_max + 1))
            plt.minorticks_on()
        # Add threshold lines and limits
        ax_p = plt.axhline(y=th_p, linestyle='--', color=COLOR_TH, linewidth=1.5)
        ax_ra = plt.axvline(x=th_ratio, linestyle='--', color=COLOR_TH, linewidth=1.5)
        ax_rb = plt.axvline(x=-th_ratio, linestyle='--', color=COLOR_TH, linewidth=1.5)
        plt.xlim(x_min, x_max)
        plt.ylim(0, y_max)
        if "_" in title:
            title = title.replace("_", " ")
        plt.title(title, fontweight="bold")
        # Set labels using iterative optimization to avoid overlaps
        objects = [ax_p, ax_ra, ax_rb]
        _set_labels(labels,
                    objects=objects,
                    th_filter=th_filter,
                    fig_format=fig_format,
                    force_points=force[0],
                    force_text=force[1],
                    force_objects=force[2],
                    box=box,
                    verbose=verbose,
                    precision=precision,
                    label_size=label_size,
                    label_bold=label_bold)
        # Add legend
        list_col = [COLOR_UP, COLOR_DOWN, COLOR_NOT_SIG]
        list_legend = ["Up", "Down", "Not Sig"]
        patches = [mpatches.Patch(color=color, label=label, hatch="o") for color, label in zip(list_col, list_legend)]
        plt.legend(handles=patches, loc=loc_legnd, frameon=False)
        # Save or show plot
        sns.despine()
        plt.tight_layout()
        title_out = title.replace(" ", "_")
        plt.xlabel(col_ratio, weight="bold")
        plt.ylabel(col_pval, weight="bold")
        ax = plt.gca()
        return ax

    @staticmethod
    def volcano_plot_ia(df_ratio_pval=None, th_filter=(0.05, 0.5), title=None,
                        col_ratio=None, col_pval=None):
        """Interactive volcano plot"""
        if title is None:
            title = "Volcano Plot for KO vs WT"
        # Adjust threshold
        th_p, th_ratio = th_filter
        if th_p < 0.5:
            th_p = -np.log10(th_p)
        colors = _color_filter(df=df_ratio_pval,
                               col_ratio=col_ratio,
                               col_pval=col_pval,
                               th_p=th_p,
                               th_ratio=th_ratio)
        dict_color_label = {COLOR_UP: "Up", COLOR_DOWN: "Down", COLOR_NOT_SIG: "Not Sig"}
        colors = [dict_color_label[color] for color in colors]
        fig = px.scatter(df_ratio_pval, hover_name="Gene_Name",
                         labels={col_ratio: col_ratio,
                                 col_pval: col_pval,
                                 "color": "status"},
                         x=col_ratio, y=col_pval,
                         color=colors, template="plotly_white", marginal_y="violin")
        fig.update_layout(title_text=title, title_x=0.5)
        fig.show()
