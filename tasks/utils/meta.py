"""
meta.yaml file module for submissions
"""

import yaml


def write_meta(data, file_path):
    with open(file_path, "w") as f:
        yaml.dump(data, f)
