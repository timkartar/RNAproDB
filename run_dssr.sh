#/bin/bash
split -l 100 <(ls ../data/cifs/) ./txts/cifs- --numeric-suffixes
for f in ./txts/*;
do
    python run_dssr.py $f &
done
