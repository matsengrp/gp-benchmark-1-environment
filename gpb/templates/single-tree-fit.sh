set -eux

trap "rm -f mmap.dat" EXIT

declare -i count=0
file=rerooted-topologies.noburnin.withbranchlengths.every100.nwk
while IFS= read -r <&3 line || [[ -n "$line" ]]; do
	tree_out=gp.$count
	TMPFILE=$(mktemp /tmp/single_tree-XXXXX.nwk)
	echo "$line" >> $TMPFILE
	gpb fit --config {{config_path}} $TMPFILE {{output_prefix}}.fasta $tree_out
	gpb pcspsurface	--config {{config_path}} $TMPFILE {{output_prefix}}.fasta $tree_out
	rm $TMPFILE
	count+=1
	if [[ $count -eq 1 ]]; then
		break
	fi
done 3< "$file"
