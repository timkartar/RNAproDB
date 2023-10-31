#!/bin/bash
sed 's/,/\n/g' ../data/nakb_prna_ids.txt > ../data/nakb_prna_ids_newlined.txt

c=1
while read p; do
        python rna_vis.py $p
        ((c++));
        echo $c
done < ../data/nakb_prna_ids_newlined.txt
