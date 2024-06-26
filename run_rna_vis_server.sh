#!/bin/bash

# PATH=$PATH:/home/aricohen/Desktop/DeepPBS/dependencies/bin
PATH=$PATH:/srv/www/deeppbs.usc.edu/deeppbs-webserver/deeppbs/dependencies/bin
PATH=$PATH:/srv/www/rnaprodb/rnaprodb_dev/external
cd /srv/www/rnaprodb/rnaprodb_dev
/home/aricohen/.conda/envs/rnascape/bin/python ./rna_vis.py $1 $2 $3 >output.log 2>error.log
