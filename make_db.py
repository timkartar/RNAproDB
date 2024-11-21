import sqlite3
import json
import os, sys
from sqlite3 import Error
import pypdb
import json
from tqdm import tqdm

def create_connection(db_file):
    """ create a database connection to a SQLite database """
    conn = None
    conn = sqlite3.connect(db_file)
    print(sqlite3.version)
    #except Error as e:
    #print(e)
    if conn:
        conn.close()
        print("closed connection")

def create_tables(f):
    conn = sqlite3.connect(f)
    cursor = conn.cursor()
    cursor.execute("drop table if exists Structures")
    sql = '''
    create table Structures(
    PDB_ID char(4) NOT NULL, 
    AUTHORS text,
    TITLE text,
    YEAR INT,
    DOI text,
    PubMed INT,
    IS_RNA_PROTEIN bool
    )
    '''

    cursor.execute(sql)
    conn.commit()
    conn.close()

def add_data(conn, table_name, data):
    query = ("insert into {} values {}".format(table_name, tuple(data)))
    print(query)
    conn.execute(query)
    conn.commit()


if __name__ == '__main__':
    dbf = "./sqldb/test.db"
    create_tables(dbf)
    conn = sqlite3.connect(dbf)  
    for item in tqdm(os.listdir("output")):
        if "pca_graph" not in item:
            continue
        else:
            try:
                prefix = item.split("_")[0]
                info = pypdb.get_info(prefix)
                data = [prefix]
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
                try:   
                    protein_count = info['rcsb_entry_info']['polymer_entity_count_protein']
                    nucleic_acid_count = info['rcsb_entry_info']['polymer_entity_count_RNA']

                    if protein_count >= 1 and nucleic_acid_count >= 1: 
                        data.append(True)
                    else:
                        data.append(False)
                except Exception as e:
                    print(e)
                    data.append("NULL")
                add_data(conn, "Structures", data)
            except Exception as e:
                print(e)
    
    conn.close()