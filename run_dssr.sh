#/bin/bash
split -l 100 <(ls ../data/cifs/) ./txts/cifs- --numeric-suffixes
for FILE in /home/raktim/rnaprodb/rnaprodb/txts/$1*;
do  
    f="$(basename -- $FILE)"
    #mkdir ./run/${f}_run
    cp *.py ./run/${f}_run/
    #cp -r _data ./run/${f}_run/
    cd ./run/${f}_run
    python run_dssr.py $FILE &
    cd /home/raktim/rnaprodb/rnaprodb/
done
