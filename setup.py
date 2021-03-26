"""Our setup script."""

import glob
from setuptools import setup

setup(
    name="gpb",
    description="Benchmarking generalized pruning",
    packages=["gpb"],
    package_data={"gpb": ["data/*"]},
    scripts=glob.glob("gpb/scripts/*.sh"),
    entry_points={"console_scripts": ["gpb=gpb.cli:safe_cli"]},
    install_requires=[
        "click_config_file",
        "csvkit",
        "jinja2",
        "plotnine",
        "seaborn",
        "seqmagick",
    ],
)
