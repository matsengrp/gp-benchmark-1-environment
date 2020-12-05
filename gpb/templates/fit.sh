set -eux

trap "rm -f mmap.dat" EXIT
gpb fit --config {{config_path}} rerooted-topologies.noburnin.nwk {{output_prefix}}.fasta gp-sbn-parameters.csv

# Train via SBN-SA and put in sa-sbn-params.csv
# python train-sbn-sa.py
#
# python plot.py

