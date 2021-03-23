# Draw a tree uniformly from the input file, and draw a uniform NNI neighbor from that file.

in_file=$1

sort -R $in_file | head -n 1 \
    | spr_neighbors --nni \
    | sort -R | head -n 1 \
    | nw_order -
