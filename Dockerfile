FROM quay.io/matsengrp/bito

RUN apt-get --allow-releaseinfo-change update
RUN apt-get install -y mrbayes libboost-graph-dev

# Make RUN commands use the bito environment:
SHELL ["/opt/conda/bin/conda", "run", "-n", "bito", "/bin/bash", "-c"]

RUN conda install -c bioconda newick_utils
RUN conda install -c conda-forge parallel

RUN git clone --recurse-submodules https://github.com/phylovi/bito.git
WORKDIR /bito

COPY . /gp-benchmark-1-environment
WORKDIR /gp-benchmark-1-environment
RUN git submodule update --init --recursive
RUN make -C spr_neighbors
RUN cp spr_neighbors/spr_neighbors $CONDA_PREFIX/bin
RUN pip install -e .
