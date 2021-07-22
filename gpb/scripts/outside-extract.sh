# An accessory script called by `outside.sh` to call outsideprob.

set -eu

outside=$1
base=$(basename $outside .sbn.csv)

gpb outsideprob _ignore_outside/rerooted-topologies.noburnin.nonsingletons.gp.sbn.csv $outside _ignore_outside/${base}.outside.csv
