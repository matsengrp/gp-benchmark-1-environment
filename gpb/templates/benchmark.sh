set -eux

trap "rm -f mmap.dat" EXIT

TREES=rerooted-topologies.noburnin.withbranchlengths.nwk

gpb benchmark --config {{config_path}} $TREES {{output_prefix}}.fasta {{output_prefix}}.bench
