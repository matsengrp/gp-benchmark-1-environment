set -eux

trap "rm -f mmap.dat" EXIT

TREES=rerooted-topologies.noburnin.nonsingletons.nwk

gpb fit --config {{config_path}} $TREES {{output_prefix}}.fasta gp-sbn-parameters.csv
gpb sa $TREES sa-sbn-parameters.csv
gpb compare gp-sbn-parameters.csv sa-sbn-parameters.csv comparison
