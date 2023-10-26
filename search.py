import pypdb

def query_by_term(term: str) -> list:
    query_results = pypdb.Query(term).search()
    with open('hirad/nakb_prna_ids.txt', 'r') as f:
        nakb_ids = f.read().split(',')
    output_list = [entry.upper() for entry in nakb_ids if entry.upper() in set(query_results)]    
    return output_list

def query_by_uniprot(uniprot: str) -> list:
   query_results = pypdb.Query(uniprot, query_type="uniprot").search()
   with open('hirad/nakb_prna_ids.txt', 'r') as f:
    nakb_ids = f.read().split(',')
   output_list = [entry.upper() for entry in nakb_ids if entry.upper() in set(query_results)]    
   return output_list

def query_by_pfam(pfam: str) -> list:
   query_results = pypdb.Query(pfam, query_type="pfam").search()
   with open('hirad/nakb_prna_ids.txt', 'r') as f:
    nakb_ids = f.read().split(',')
   output_list = [entry.upper() for entry in nakb_ids if entry.upper() in set(query_results)]    
   return output_list

def query_by_pmid(pmid: int) -> list:
   query_results = pypdb.Query(pmid, query_type="PubmedIdQuery").search()
   with open('hirad/nakb_prna_ids.txt', 'r') as f:
    nakb_ids = f.read().split(',')
   output_list = [entry.upper() for entry in nakb_ids if entry.upper() in set(query_results)]    
   return output_list


 
if __name__ == "__main__":
    print(query_by_pmid(27499440))
