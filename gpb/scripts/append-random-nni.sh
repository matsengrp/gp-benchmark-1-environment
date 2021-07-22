# Append a random NNI to a file, making sure that it's not in the file already.
# If the file has every NNI for every tree, this script will loop infinitely.
# Name the file with the NNI according to the SHA of the new tree.

set -eu

in_file=$1
out_prefix=$2

while true; do
    neighbor_tree=$(random-nni.sh $in_file)
    grep -q "$neighbor_tree" $in_file || break
done

out_file=${out_prefix}$(echo -n "$neighbor_tree" | sha256sum | cut -f 1 -d ' ').nwk

cp $in_file $out_file
echo $neighbor_tree >> $out_file

