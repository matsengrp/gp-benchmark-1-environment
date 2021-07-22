# A pipeline for learning about the GP probability for GPCSPs outside of the support.

set -eu

fasta_path=$1
target_nni_count=$2

rm -rf _ignore_outside
mkdir -p _ignore_outside

# ### Start by generating trees ###
# We do a little work here because we will often generate the same NNI multiple times.
# Because we're using the SHA of the generated tree as a filename, the same NNI will result in the same filename.
generate_count=$(($target_nni_count + $target_nni_count/5))

seq $generate_count | parallel append-random-nni.sh rerooted-topologies.noburnin.nonsingletons.nwk "_ignore_outside/"

actual_generated_count=$(ls _ignore_outside/*nwk | wc -l)
cull_count=$(($actual_generated_count - $target_nni_count))
rm $(ls _ignore_outside/*nwk | tail -n $cull_count)

test $target_nni_count -eq $(ls _ignore_outside/*nwk | wc -l) || {
    echo "Didn't get the right number of NNIs."
    exit 1
}

echo "Successfully generated $target_nni_count NNIs."


# ### Perform GP
ls _ignore_outside/*nwk | parallel outside-fit.sh $fasta_path
outside-fit.sh $fasta_path rerooted-topologies.noburnin.nonsingletons.nwk

echo "Successfully fit GP."


# ### Gather the results ###
ls _ignore_outside/*.sbn.csv | grep -v rooted | parallel outside-extract.sh
csvstack _ignore_outside/*outside.csv > all_outside.csv
gpb addmeta _ignore_outside/rerooted-topologies.noburnin.nonsingletons.gp.sbn.csv _ignore_outside/rerooted-topologies.noburnin.nonsingletons.gp.sbn.meta.csv
for out_path in outside.svg outside.pdf; do
    gpb outsideplot all_outside.csv _ignore_outside/rerooted-topologies.noburnin.nonsingletons.gp.sbn.meta.csv $out_path
done
