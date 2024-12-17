#!/bin/bash

PATH=$PATH:/srv/www/deeppbs.usc.edu/deeppbs-webserver/deeppbs/dependencies/bin
PATH=$PATH:/srv/www/rnaprodb/rnaprodb_dev/external

cd /srv/www/rnaprodb/rnaprodb_dev

source /opt/anaconda3/etc/profile.d/conda.sh

conda activate rnascape
nohup python auto_update.py
