set -eux

trap "rm -f mmap.dat" EXIT

UNIFTREES=rerooted-topologies.unif.noburnin.withbranchlengths.nwk
EXPTREES=rerooted-topologies.exp.noburnin.withbranchlengths.nwk
UNIQTREES=rerooted-topologies.unif.noburnin.uniq.nwk

for i in {1..{{benchmark_iters}}}
do
    gpb timingbenchmark --config {{config_path}} $UNIFTREES {{output_prefix}}.fasta {{output_prefix}}.bench.timing
done

for i in {1..{{benchmark_iters}}}
do
    gpb timingbenchmark --config {{config_path}} $UNIQTREES {{output_prefix}}.fasta {{output_prefix}}.bench.timing.uniq
done

gpb estimationbenchmark --config {{config_path}} $UNIFTREES $EXPTREES {{output_prefix}}.fasta {{output_prefix}}.bench.estimation

