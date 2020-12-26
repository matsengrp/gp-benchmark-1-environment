set -eux

trap "rm -f mmap.dat" EXIT

TREES=rerooted-topologies.noburnin.nonsingletons.nwk

gpb fit --config {{config_path}} $TREES {{output_prefix}}.fasta gp
gpb sa $TREES sa-sbn-parameters
gpb compare gp.sbn.csv sa-sbn-parameters.csv sa-sbn-parameters.subsplit.csv {{output_prefix}}.comparison
