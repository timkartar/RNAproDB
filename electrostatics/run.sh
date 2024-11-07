#!/bin/bash
export PATH="/srv/www/rnaprodb/rnaprodb_dev/electrostatics/TABI-PB/build/bin:/srv/www/rnaprodb/rnaprodb_dev/electrostatics/pnabind/share:/home/aricohen/.conda/envs/rnascape/bin/pdb2pqr:$PATH"
/home/aricohen/.conda/envs/rnascape/bin/pdb2pqr --ff=AMBER ./$1.pdb $1.pqr
/home/aricohen/.conda/envs/rnascape/bin/python ./pnabind/pnabind/structure/run_tabipb.py $1
/home/aricohen/.conda/envs/rnascape/bin/python vis.py $1
