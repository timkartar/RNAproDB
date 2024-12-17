#!/bin/bash

# PATH=$PATH:/home/aricohen/Desktop/DeepPBS/dependencies/bin
PATH=$PATH:/srv/www/rnaprodb.usc.edu/DeepPBS/dependencies/bin
PATH=$PATH:/srv/www/rnaprodb.usc.edu/rnaprodb/rnaprodb_dev/external

cd /srv/www/rnaprodb.usc.edu/rnaprodb/rnaprodb_dev
rm -f output.log
rm -f $1-assembly1.pdb
rm -f $1-assembly1.hb2
rm -f error.log
/srv/www/rnaprodb.usc.edu/conda/envs/rnaprodb/bin/python ./rna_vis.py $1 $2 $3
# >output.log 2>error.log
cd -
