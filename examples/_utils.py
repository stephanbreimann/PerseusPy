"""
This is a script with utility functions
"""
import os
import platform


def _folder_path(super_folder, folder_name):
    """Modification of separator (OS depending)"""
    path = os.path.join(super_folder, folder_name + SEP)
    return path


# Folder
SEP = "\\" if platform.system() == "Windows" else "/"

FOLDER_PROJECT = f"{SEP}".join(os.path.dirname(os.path.abspath(__file__)).split("/")[0:-1]) + SEP
FOLDER_DATA = _folder_path(FOLDER_PROJECT, 'data')
