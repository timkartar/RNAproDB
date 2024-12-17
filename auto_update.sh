#!/bin/bash

PATH=$PATH:/srv/www/rnaprodb.usc.edu/DeepPBS/dependencies/bin
PATH=$PATH:/srv/www/rnaprodb.usc.edu/rnaprodb/rnaprodb_dev/external


cd /srv/www/rnaprodb.usc.edu/rnaprodb/rnaprodb_dev

source /srv/www/rnaprodb.usc.edu/conda/etc/profile.d/conda.sh

conda activate rnaprodb
nohup python auto_update.py
