set -eu

# Running iqtree and rooting ML tree on {{reroot_number}}, which is {{reroot_name}}
iqtree2 -s {{output_prefix}}.fasta -m JC69 -pre {{output_prefix}}.ml -redo -quiet
cat {{output_prefix}}.ml.treefile | nw_reroot - {{reroot_number}} > {{output_prefix}}.ml.nwk

seqmagick convert {{nexus_alignment_path}} {{output_prefix}}.fasta
