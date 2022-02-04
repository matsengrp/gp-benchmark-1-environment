set -eux

trap "rm -f mmap.dat" EXIT

TREES=rerooted-topologies.noburnin.nonsingletons.nwk

gpb fit --config {{config_path}} $TREES {{output_prefix}}.fasta gp
#gpb pcspsurface --config {{config_path}} $TREES {{output_prefix}}.fasta gp
gpb sa $TREES sa.sbn
gpb compare gp.sbn.csv sa.sbn.csv sa.sbn.subsplit.csv {{output_prefix}}.comparison

for out_path in {{output_prefix}}.perpcsp.plot.svg {{output_prefix}}.perpcsp.plot.pdf; do
	gpb pcspoptplot gp.perpcsp_llh_from_opt.csv $out_path
done

# reading in path for likelihood surfaces obtained after optimization
#for out_path in {{output_prefix}}.perpcsp.surface.svg {{output_prefix}}.perpcsp.surface.pdf; do
#	gpb pcspsurfaceplot gp.perpcsp_llh_surface.csv $out_path
#done
