set -eu

# rerooting on {{reroot_number}}, which is {{reroot_name}}
awk '$1~/tree/ {print $NF}' {{output_prefix}}.t | nw_reroot - {{reroot_number}} | nw_order - | tail -n +{{burnin_samples}} > rerooted-topologies.noburnin.withbranchlengths.nwk

NUM=100
awk -v NUM=$NUM 'NR % NUM == 0' rerooted-topologies.noburnin.withbranchlengths.nwk | sort -nr > rerooted-topologies.noburnin.withbranchlengths.every100.nwk 
#sort rerooted-topologies.noburnin.withbranchlengths.nwk | uniq -c | sort -nr | awk '$1 > 1 {print $2}' > rerooted-topologies.noburnin.withbranchlengths.nonsingletons.nwk

seqmagick convert {{nexus_alignment_path}} {{output_prefix}}.fasta
