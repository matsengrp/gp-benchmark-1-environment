set -eu

# rerooting on {{reroot_number}}, which is {{reroot_name}}
awk '$1~/tree/ {print $NF}' {{output_prefix}}.t | nw_topology - | nw_reroot - {{reroot_number}} | nw_order - \
    | uniq -c \
    | sed -e "s/^[ ]*//" -e "s/[ ]/\t/" \
    > rerooted-topology-sequence.tab
