"""Our setup script."""

import glob
from setuptools import setup

setup(
    name="gpb",
    description="Benchmarking generalized pruning",
    packages=["gpb"],
    package_data={"gpb": ["data/*"]},
    scripts=glob.glob("gpb/scripts/*.sh"),
    entry_points={"console_scripts": ["gpb=gpb.cli:cli"]},
    install_requires=[
        "click-config-file",
        "csvkit",
        "dendropy",
        "jinja2",
        "plotnine==0.6.0",
        "pylint",
        "pytest",
    ],
)
