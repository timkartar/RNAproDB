#!/bin/bash
while read p; do
    python run.py $p
    ./run.sh full_$p
    #./run.sh na_$p
    #./run.sh pro_$p
    rm pro_*
    rm na_*
    rm full_*
done <$1
