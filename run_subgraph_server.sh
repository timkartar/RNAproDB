#!/bin/bash
PATH=$PATH:/srv/www/deeppbs.usc.edu/deeppbs-webserver/deeppbs/dependencies/bin
cd /srv/www/rnaprodb/rnaprodb_dev
/home/aricohen/.conda/envs/rnascape/bin/python ./get_subgraph.py $1 $2 $3
# >output.log 2>error.log
