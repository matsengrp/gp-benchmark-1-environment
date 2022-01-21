set -eux

trap "rm -f mmap.dat" EXIT

TREEFILE={{output_prefix}}.ml.nwk
tree_out=gp.ml
gpb fit --config {{config_path}} $TREEFILE {{output_prefix}}.fasta $tree_out
gpb pcspsurface --config {{config_path}} $TREEFILE {{output_prefix}}.fasta $tree_out
