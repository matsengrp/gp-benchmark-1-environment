# gp-benchmark-1-environment

[![Docker Repository on Quay](https://quay.io/repository/matsengrp/gp-benchmark-1-environment/status "Docker Repository on Quay")](https://quay.io/repository/matsengrp/gp-benchmark-1-environment)

## Dependencies

We assume that this will be run in the [libsbn](https://github.com/phylovi/libsbn) conda environment.

Python dependencies will be installed via pip with the command below.
We also require:

* [MrBayes](https://nbisweden.github.io/MrBayes/index.html)
* [Newick utilities](http://cegg.unige.ch/newick_utils)
* [spr_neighbors](https://github.com/cwhidden/spr_neighbors)
* [GNU Parallel](https://www.gnu.org/software/parallel/)

Please see Dockerfile for complete install details.


## Install

    pip install -e .
