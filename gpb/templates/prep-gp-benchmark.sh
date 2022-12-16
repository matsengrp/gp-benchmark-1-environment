set -eu

# rerooting on {{reroot_number}}, which is {{reroot_name}}
# doing this on the from the MB runs with a uniform prior to load into GP
# and also on the trees from MB  with an exponential prior to assess VBPI fit
UNIFPATH=../MrBayes/unif
EXPPATH=../MrBayes/exp
awk '$1~/tree/ {print $NF}' $UNIFPATH/{{output_prefix}}.t | nw_reroot - {{reroot_number}} | nw_order - | tail -n +{{burnin_samples}} > rerooted-topologies.unif.noburnin.withbranchlengths.nwk

awk '$1~/tree/ {print $NF}' $EXPPATH/{{output_prefix}}.t | nw_reroot - {{reroot_number}} | nw_order - | tail -n +{{burnin_samples}} > rerooted-topologies.exp.noburnin.withbranchlengths.nwk

# outputting unique topologies to assess how much faster it is build engine
nw_topology rerooted-topologies.unif.noburnin.withbranchlengths.nwk | sort | uniq > rerooted-topologies.unif.noburnin.uniq.nwk

seqmagick convert {{nexus_alignment_path}} {{output_prefix}}.fasta
