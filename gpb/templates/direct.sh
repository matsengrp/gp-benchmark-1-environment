set -eu

trap "rm -f mmap.dat" EXIT

mkdir -p _ignore/trees
gpb treeexport --config config.json rerooted-topologies.noburnin.nonsingletons.nwk {{output_prefix}}.fasta gp.sbn.csv

for i in _ignore/trees/*.nwk; do
    echo $i | grep -q deroot && continue
    if test -s $i; then
        nw_reroot -d $i > $(dirname $i)/$(basename $i .nwk).deroot.nwk
    fi
done

gpb treemarginal "_ignore/trees/*.deroot.nwk" {{output_prefix}}.fasta pcsp-marginals.csv
gpb comparedirect pcsp-marginals.csv gp.prior.csv sa.sbn.csv direct.sbn
