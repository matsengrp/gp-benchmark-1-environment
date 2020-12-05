set -eux

trap "rm -f mmap.dat" EXIT

gpb fit --config {{config_path}} rerooted-topologies.noburnin.nwk {{output_prefix}}.fasta gp-sbn-parameters.csv
gpb sa rerooted-topologies.noburnin.nwk sa-sbn-parameters.csv
gpb compare gp-sbn-parameters.csv sa-sbn-parameters.csv comparison
