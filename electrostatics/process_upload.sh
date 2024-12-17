#!/bin/bash
touch /srv/www/rnaprodb.usc.edu/rnaprodb/rnaprodb_dev/output/electrostatics/ip_$1.txt
chmod 777 /srv/www/rnaprodb.usc.edu/rnaprodb/rnaprodb_dev/output/electrostatics/ip_$1.txt
/srv/www/rnaprodb.usc.edu/conda/envs/rnaprodb/bin/python run.py $1
./run.sh full_$1
./run.sh na_$1
./run.sh pro_$1
chmod 777 /srv/www/rnaprodb.usc.edu/rnaprodb/rnaprodb_dev/output/electrostatics/full_$1.ply
chmod 777 /srv/www/rnaprodb.usc.edu/rnaprodb/rnaprodb_dev/output/electrostatics/pro_$1.ply
chmod 777 /srv/www/rnaprodb.usc.edu/rnaprodb/rnaprodb_dev/output/electrostatics/na_$1.ply
rm pro_*
rm na_*
rm full_*
rm /srv/www/rnaprodb.usc.edu/rnaprodb/rnaprodb_dev/output/electrostatics/ip_$1.txt
