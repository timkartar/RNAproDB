import sqlite3
import os, django


def import_SQL():
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "backend.settings")
    django.setup()
    from search.models import pypdbObject
    conn = sqlite3.connect('../sqldb/test.db')
    c = conn.cursor()
    c.execute("SELECT * FROM Structures")
    rows = c.fetchall()
    pypdbObject.objects.all().delete()
    for row in rows:
        # Convert 'NULL' string to None for nullable fields
        authors = None if row[1] == 'NULL' else row[1]
        title = None if row[2] == 'NULL' else row[2]
        year = None if row[3] == 'NULL' else row[3]
        doi = None if row[5] == 'NULL' else row[5]
        pubmed = None if row[4] == 'NULL' else row[4]
        is_rna_protein = None if row[6] == 'NULL' else row[6]

        pypdbObject.objects.create(
            id=row[0], 
            authors=authors, 
            title=title, 
            year=year, 
            pubmed=pubmed, 
            doi=doi,
            is_rna_protein=is_rna_protein
        )
    conn.commit()
    conn.close()

def main():
    import_SQL()

if __name__ == '__main__':
    main()