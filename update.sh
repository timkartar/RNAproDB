#!/bin/bash
rm new_ids.txt
python download_ids.py > new_ids.txt
chmod 777 new_ids.txt

echo "All considered ids":
cat ../rnaprodb_frontend/build/ids.txt | sed "s/,/\n/g" | wc -l
echo "All fetched ids":
cat new_ids.txt | wc -l
for id in $(grep -v -F -x -f <(cat ../rnaprodb_frontend/build/ids.txt | sed "s/,/\n/g") new_ids.txt); do
	./run_rna_vis_server.sh $id
done

rm ../rnaprodb_frontend/build/ids_in_database.txt
ls ./output/*_pca_* | cut -d \/ -f 3 | cut -d _ -f 1  | awk '$1=$1' RS= OFS=, | cat > ../rnaprodb_frontend/build/ids_in_database.txt
chmod 777 ../rnaprodb_frontend/build/ids_in_database.txt

mv ../rnaprodb_frontend/build/ids.txt ../rnaprodb_frontend/build/ids_last_update.txt
cat new_ids.txt | awk '$1=$1' RS= OFS=, > ../rnaprodb_frontend/build/ids.txt
chmod 777 ../rnaprodb_frontend/build/ids.txt

chmod -R 777 ../rnaprodb_frontend/build/cifs/
