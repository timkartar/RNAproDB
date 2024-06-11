from urllib.request import urlretrieve
import gzip 
import shutil
import pypdb
from pypdb.clients.search.search_client import perform_search, perform_search_with_graph, ReturnType, QueryGroup, LogicalOperator
from pypdb.clients.search.operators import text_operators
import sqlite3
import subprocess
import os

def get_pdb_file(pdb_id):
    base_script_path = '/mnt/c/Users/hosse/Desktop/Github/frontend_master'
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
    with open('/mnt/c/Users/hosse/Desktop/Github/frontend_master/public/ids.txt', 'r+') as f:
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
    dbf = "./sqldb/test.db"
    conn = sqlite3.connect(dbf)
    info = pypdb.get_info(pdb_id)
    data = [pdb_id]
    data.append(",".join(info['citation'][0]['rcsb_authors']))
    data.append(info['citation'][0]['title'].replace("\'", "single_quote"))
    try:
        data.append(info['citation'][0]['year'])
    except:
        data.append("NULL")
    try:
        data.append(info['citation'][0]['pdbx_database_id_pub_med'])
    except:
        data.append("NULL")
    try:
        data.append(info['citation'][0]['pdbx_database_id_doi'])
    except:
        data.append("NULL")
    query = ("insert into {} values {}".format("Structures", tuple(data)))
    conn.execute(query)
    conn.commit()
    conn.close()
    return 'Django database updated'

def run_django_migrations():
    try:
        manage_py_path = '/mnt/c/Users/hosse/Desktop/Github/backend_master/backend/manage.py'
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
    subprocess.run(['python', '/mnt/c/Users/hosse/Desktop/Github/backend_master/rna_vis.py', pdb_id], check=True)
    return 'new PDB has been added to database'

if __name__ == "__main__":
    pdb_id = '1fjg'
    print(check_pdb(pdb_id))

    




