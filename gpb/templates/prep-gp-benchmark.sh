set -eu

# rerooting on {{reroot_number}}, which is {{reroot_name}}
awk '$1~/tree/ {print $NF}' {{output_prefix}}.t | nw_reroot - {{reroot_number}} | nw_order - | tail -n +{{burnin_samples}} > rerooted-topologies.noburnin.withbranchlengths.nwk

nw_topology rerooted-topologies.noburnin.withbranchlengths.nwk | uniq -c | sort -nr > rerooted-topologies.noburnin.nwk

awk '$1 > 1 {print $2}' rerooted-topologies.noburnin.nwk > rerooted-topologies-noburnin-nonsingletons.nwk

seqmagick convert {{nexus_alignment_path}} {{output_prefix}}.fasta
