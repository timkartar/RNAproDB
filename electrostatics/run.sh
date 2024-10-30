#!/bin/bash
export PATH="/srv/www/rnaprodb/rnaprodb_dev/electrostatics/TABI-PB/build/bin:/srv/www/rnaprodb/rnaprodb_dev/electrostatics/pnabind/share:$PATH"
pdb2pqr --ff=AMBER ./$1.pdb $1.pqr
python ./pnabind/pnabind/structure/run_tabipb.py $1
# python vis.py $1
