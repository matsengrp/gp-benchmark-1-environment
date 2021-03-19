set -eux

trap "rm -f mmap.dat" EXIT

TREES=rerooted-topologies.noburnin.nonsingletons.nwk

gpb fit --config {{config_path}} $TREES {{output_prefix}}.fasta gp
gpb sa $TREES sa.sbn
gpb compare gp.sbn.csv sa.sbn.csv sa.sbn.subsplit.csv {{output_prefix}}.comparison
