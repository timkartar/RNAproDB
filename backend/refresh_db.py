from urllib.request import urlretrieve
import gzip 
import shutil
import pypdb
from pypdb.clients.search.search_client import perform_search, perform_search_with_graph, ReturnType, QueryGroup, LogicalOperator
from pypdb.clients.search.operators import text_operators
import sqlite3
import subprocess
import os
import django

backend =  os.path.dirname(os.path.abspath(__file__))
frontend = "/srv/www/rnaprodb/rnaprodb_frontend/"

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "backend.settings")
django.setup()
from search.models import pypdbObject


def get_pdb_file(pdb_id):
    base_script_path = frontend
    pdb_id = pdb_id.strip().upper()
    if not pdb_id.isalnum():
        return 'error', 'Bad PDB ID'
    out_gz_path = '{}/public/cifs/{}-assembly1.cif.gz'.format(base_script_path, pdb_id)
    out_cif_path = '{}/public/cifs/{}-assembly1.cif'.format(base_script_path, pdb_id)
    url = 'https://files.rcsb.org/download/{}-assembly1.cif.gz'.format(pdb_id)
    urlretrieve(url, out_gz_path)
    with gzip.open(out_gz_path, 'rb') as f_in:
        with open(out_cif_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    filepath = f'cifs/{pdb_id}-assembly1.cif'
    # if not os.path.exists(filepath): # if path doesn't exist
    #     return 'error', 'Bad PDB ID'
        # while not is_download_complete():
    #     time.sleep(1)
    return 'PDB file downloaded successfully', filepath

def process_pdb(pbd_id):
    subprocess.call(['python', 'run_dssr.py', f'{pbd_id}'])
    with open('nakb_prna_ids.txt', 'r+') as f:
        lines = f.readlines()
        text_to_append = pbd_id + ','
        if lines:
            lines[-1] = lines[-1].rstrip('\n')  
            lines[-1] += text_to_append
        else:
            lines.append(text_to_append) 
        f.seek(0)
        f.writelines(lines)
        f.truncate()  
    with open(frontend + '/public/ids.txt', 'r+') as f:
        lines = f.readlines()
        text_to_append = pbd_id + ','
        if lines:
            lines[-1] = lines[-1].rstrip('\n')  
            lines[-1] += text_to_append
        else:
            lines.append(text_to_append) 
        f.seek(0)
        f.writelines(lines)
        f.truncate()  
    return 'DSSR ran successfully'

def populate_db(pdb_id):
    dbf = "../sqldb/test.db"
    conn = sqlite3.connect(dbf)
    info = pypdb.get_info(pdb_id)
    data = [pdb_id]
    
    authors = ",".join(info['citation'][1]['rcsb_authors'])
    title = info['citation'][0]['title'].replace("\'", "single_quote")
    year = info['citation'][0]['year']
    pubmed = info['citation'][0]['pdbx_database_id_pub_med']
    doi = info['citation'][0]['pdbx_database_id_doi']

    try:
        data.append(authors)
    except:
        data.append("NULL")
    try:
        data.append(title)
    except:
        data.append("NULL")
    try:
        data.append(year)
    except:
        data.append("NULL")
    try:
        data.append(pubmed)
    except:
        data.append("NULL")
    try:
        data.append(doi)
    except:
        data.append("NULL")

    query = ("insert into {} values {}".format("Structures", tuple(data)))
    conn.execute(query)
    conn.commit()
    conn.close()
    pypdbObject.objects.create(
            id=pdb_id,
            authors=authors,
            title=title,
            year=year,
            pubmed=pubmed,
            doi=doi
        )
    return 'Django database updated'

def run_django_migrations():
    try:
        manage_py_path = backend + '/manage.py'
        subprocess.run(['python', manage_py_path, 'makemigrations'], check=True)

        subprocess.run(['python', manage_py_path, 'migrate'], check=True)

        print("Migrations were successfully applied.")

    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running migrations: {e}")

def check_pdb(pdb_id):
    pdb_id = pdb_id.strip().upper()
    if not pdb_id.isalnum():
        return 'Bad PDB ID'
    queries = []
    # number of nucleic acid entities greater than zero
    queries.append(pypdb.text_operators.ComparisonOperator(value = 0, attribute = "rcsb_assembly_info.polymer_entity_count_RNA", comparison_type = pypdb.text_operators.ComparisonType.GREATER))
    queries.append(pypdb.text_operators.ComparisonOperator(value = 0, attribute = "rcsb_assembly_info.polymer_entity_count_protein", comparison_type = pypdb.text_operators.ComparisonType.GREATER))
    results = perform_search_with_graph(
    query_object=QueryGroup(
        logical_operator=LogicalOperator.AND,
        queries=queries
    ),
    return_type=ReturnType.ENTRY
    )
    if pdb_id not in results:
        return 'This PDB ID does not contain RNA and protein'
    get_pdb_file(pdb_id)
    process_pdb(pdb_id)
    populate_db(pdb_id)
    run_django_migrations()
    subprocess.run(['python', backend + '/rna_vis.py', pdb_id], check=True)
    return 'new PDB has been added to database'

if __name__ == "__main__":
    import sys
    pdb_id = sys.argv[1]
    #  print(check_pdb(pdb_id))
    populate_db(pdb_id)
    run_django_migrations()
    




