#!/bin/bash
touch /home/db/rnaprodb/output/electrostatics/ip_$1.txt
chmod 777 /home/db/rnaprodb/output/electrostatics/ip_$1.txt
/home/aricohen/.conda/envs/rnascape/bin/python run.py $1
./run.sh full_$1
./run.sh na_$1
./run.sh pro_$1
chmod 777 /home/db/rnaprodb/output/electrostatics/full_$1.ply
chmod 777 /home/db/rnaprodb/output/electrostatics/pro_$1.ply
chmod 777 /home/db/rnaprodb/output/electrostatics/na_$1.ply
rm pro_*
rm na_*
rm full_*
rm /home/db/rnaprodb/output/electrostatics/ip_$1.txt
