set -eu

# rerooting on {{reroot_number}}, which is {{reroot_name}}
awk '$1~/tree/ {print $NF}' {{output_prefix}}.t | nw_topology - | nw_reroot - {{reroot_number}} | nw_order - | tail -n +{{burnin_samples}} > ds1.rerooted-topologies.noburnin.nwk
