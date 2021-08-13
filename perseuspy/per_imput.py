"""
This is a script for imputation methods in Perseus pipeline
"""
import pandas as pd

from perseuspy.per_base import PerseusBase


# I Helper Functions


# II Main Functions
class PerseusImputation(PerseusBase):
    """Class for Perseus analysis"""
    def __init__(self, **kwargs):
        PerseusBase.__init__(self, **kwargs)

    # TODO implement imputation method

