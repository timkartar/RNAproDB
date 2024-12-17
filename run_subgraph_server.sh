#!/bin/bash
PATH=$PATH:/srv/www/rnaprodb.usc.edu/DeepPBS/dependencies/bin
PATH=$PATH:/srv/www/rnaprodb.usc.edu/rnaprodb/rnaprodb_dev/external

cd /srv/www/rnaprodb.usc.edu/rnaprodb/rnaprodb_dev
/srv/www/rnaprodb.usc.edu/conda/envs/rnaprodb/bin/python ./get_subgraph.py $1 $2 $3
# >output.log 2>error.log
