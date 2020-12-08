FROM phylovi/libsbn

RUN apt-get update && apt-get install -y mrbayes

# Make RUN commands use the libsbn environment:
SHELL ["/opt/conda/bin/conda", "run", "-n", "libsbn", "/bin/bash", "-c"]

RUN conda install -c bioconda newick_utils

RUN git clone --recurse-submodules https://github.com/phylovi/libsbn.git
WORKDIR /libsbn
RUN git checkout d7a67de7322d36357f5b465c4de47eb3acb5c532
RUN make -j 20 test

COPY . /gp-benchmark-1-environment
WORKDIR /gp-benchmark-1-environment
RUN pip install -e .
