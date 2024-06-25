#!/bin/bash
PATH=$PATH:/home/aricohen/Desktop/DeepPBS/dependencies/bin
cd /home/aricohen/Desktop/django-react-rnaprodb/rnaprodb_dev
/home/aricohen/anaconda3/envs/RNAproDB/bin/python ./get_subgraph.py $1 $2 $3
# >output.log 2>error.log
