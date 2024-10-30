#!/bin/bash
while read p; do
    python run.py $p
    ./run.sh full_$p
    #./run.sh na_$p
    #./run.sh pro_$p
    #rm *$p*
done <$1
