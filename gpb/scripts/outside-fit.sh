# An accessory script called by `outside.sh` to do GP fitting.

set -eu

fasta_path=$1
tree_path=$2

base=_ignore_outside/$(basename $tree_path .nwk).gp
mmap_path=$base.mmap.dat

trap "rm -f $mmap_path" EXIT

gpb fit --tol 0.0001 --max-iter 100 --mmap-path $mmap_path $tree_path $fasta_path $base
