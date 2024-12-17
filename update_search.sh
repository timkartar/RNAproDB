#!/bin/bash

PATH=$PATH:/srv/www/deeppbs.usc.edu/deeppbs-webserver/deeppbs/dependencies/bin
PATH=$PATH:/srv/www/rnaprodb/rnaprodb_dev/external

source activate rnascape

./update_nakb_ids.py
python downloadPDBThumbnails.py
python make_db.py

cd /srv/www/rnaprodb/rnaprodb_dev/backend
python importSQL.py
