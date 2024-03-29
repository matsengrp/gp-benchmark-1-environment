set -eux

trap "rm -f mmap.dat" EXIT

TREES=rerooted-topologies.noburnin.withbranchlengths.nwk
gpb fit --config {{config_path}} $TREES {{output_prefix}}.fasta gp
gpb pcspsurface --config {{config_path}} $TREES {{output_prefix}}.fasta gp

# reading in path for likelihood surfaces obtained after optimization
for out_path in {{output_prefix}}.perpcsp.surface.svg {{output_prefix}}.perpcsp.surface.pdf; do
	gpb pcspsurfaceplot gp.perpcsp_llh_surface.csv $out_path
done
