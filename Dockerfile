FROM phylovi/libsbn

RUN apt-get update && apt-get install -y mrbayes

# Make RUN commands use the libsbn environment:
SHELL ["/opt/conda/bin/conda", "run", "-n", "libsbn", "/bin/bash", "-c"]

RUN conda install -c bioconda newick_utils

RUN git clone --recurse-submodules https://github.com/phylovi/libsbn.git
WORKDIR /libsbn
RUN git checkout 04ee68371638830e4d7302a214a72e350d3c8b63
RUN make -j 20 test

COPY . /gp-benchmark-1-environment
WORKDIR /gp-benchmark-1-environment
RUN pip install -e .
