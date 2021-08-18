"""
This is a script for imputation methods in Perseus pipeline

References
----------
[1] Liu M. and Dongre A., Proper imputation of missing values in proteomics datasets
    for differential expression analysis. Briefing in Bioinformatics (2020)
[2] Lenz M., et al. Missing values imputation in proximity extension assay-based targeted proteomics data.
    PloS One (2020)
[3] Palstrom N. B., Matthiesen R., and Beck H. C., Data Imputation in Merged Isobaric Labeling-Based
    Relative Quantification Datasets. Methods in Molecular Biology (2020)

Notes
-----
8 datasets downloaded for [1] from PRIDE :
    PXD000484
    PXD000485
    PXD004816
    PXD005065
    PXD005833
    PXD006401
    PXD006401
    PXD009933
"""
import pandas as pd

from perseuspy.per_base import PerseusBase


# I Helper Functions


# II Main Functions
class PerseusImputation(PerseusBase):
    """Class for Perseus analysis"""
    def __init__(self, **kwargs):
        PerseusBase.__init__(self, **kwargs)

    # TODO implement imputation methods

