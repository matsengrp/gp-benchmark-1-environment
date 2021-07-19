FROM quay.io/matsengrp/libsbn

RUN apt-get update && apt-get install -y mrbayes

# Make RUN commands use the libsbn environment:
SHELL ["/opt/conda/bin/conda", "run", "-n", "libsbn", "/bin/bash", "-c"]

RUN conda install -c bioconda newick_utils

RUN git clone --recurse-submodules https://github.com/phylovi/libsbn.git
WORKDIR /libsbn
RUN make -j 20 test

COPY . /gp-benchmark-1-environment
WORKDIR /gp-benchmark-1-environment
RUN git submodule update --init --recursive
RUN make -C spr_neighbors
RUN cp spr_neighbors/spr_neighbors $CONDA_PREFIX/bin
RUN pip install -e .
