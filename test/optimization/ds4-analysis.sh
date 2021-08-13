set -eu

trap "rm -f mmap.dat" EXIT

CONFIG=optim_config.json
DATA=ds4

gpb template basic.mb $CONFIG $DATA.mb
mb $DATA.mb

gpb template prep-gp.sh $CONFIG $DATA-gp-prep.sh
bash $DATA-gp-prep.sh

gpb template fit.sh $CONFIG $DATA-gp-fit.sh
bash $DATA-gp-fit.sh
