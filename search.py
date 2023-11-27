import pypdb
from pypdb.clients.search.search_client import perform_search, perform_search_with_graph, ReturnType, QueryGroup, LogicalOperator
from pypdb.clients.search.operators import text_operators
import time

# def query_by_term_control(term: str) -> list:
#     query_results = pypdb.Query(term).search()
#     with open('nakb_prna_ids.txt', 'r') as f:
#         nakb_ids = f.read().split(',')
#     output_list = [entry.upper() for entry in nakb_ids if entry.upper() in set(query_results)]    
#     return output_list

def query_by_term(term: str) -> list:
   query1 = pypdb.text_operators.ComparisonOperator(value = 0, attribute = "rcsb_assembly_info.polymer_entity_count_RNA", comparison_type = pypdb.text_operators.ComparisonType.GREATER)
   query2 = pypdb.text_operators.ComparisonOperator(value = 0, attribute = "rcsb_assembly_info.polymer_entity_count_protein", comparison_type = pypdb.text_operators.ComparisonType.GREATER)
   search_operator = text_operators.DefaultOperator(value = term)
   rna_protein_query = QueryGroup(queries = [query1, query2], logical_operator = LogicalOperator.AND)
   search_query = QueryGroup(queries = [rna_protein_query, search_operator], logical_operator = LogicalOperator.AND)
   return_type = ReturnType.ENTRY
   with open('nakb_prna_ids.txt', 'r') as f:
    nakb_ids = f.read().split(',')
   query_results = perform_search(search_query, return_type)
   output_list = [entry.upper() for entry in nakb_ids if entry.upper() in set(query_results)]
   return output_list


# def query_by_term_exp(term: str) -> list:
#    query1 = f"polymer_entity_count_rna > 0"
#    query2 = f"polymer_entity_count_protein > 0"
#    query3 = f"keywords:{term}"
#    combined_query = f"{query1} AND {query2} AND {query3}"
#    query_results = pypdb.Query(combined_query).search()
#    with open('hirad/nakb_prna_ids.txt', 'r') as f:
#        nakb_ids = f.read().split(',')
#    output_list = [entry.upper() for entry in nakb_ids if entry.upper() in set(query_results)]
#    return output_list

# def query_by_uniprot(uniprot: str) -> list:
#    query_results = pypdb.Query(uniprot, query_type="uniprot").search()
#    with open('hirad/nakb_prna_ids.txt', 'r') as f:
#     nakb_ids = f.read().split(',')
#    output_list = [entry.upper() for entry in nakb_ids if entry.upper() in set(query_results)]    
#    return output_list

# def query_by_pfam(pfam: str) -> list:
#    query_results = pypdb.Query(pfam, query_type="pfam").search()
#    with open('hirad/nakb_prna_ids.txt', 'r') as f:
#     nakb_ids = f.read().split(',')
#    output_list = [entry.upper() for entry in nakb_ids if entry.upper() in set(query_results)]    
#    return output_list

# def query_by_pmid(pmid: int) -> list:
#    query_results = pypdb.Query(pmid, query_type="PubmedIdQuery").search()
#    with open('hirad/nakb_prna_ids.txt', 'r') as f:
#     nakb_ids = f.read().split(',')
#    output_list = [entry.upper() for entry in nakb_ids if entry.upper() in set(query_results)]    
#    return output_list


if __name__ == "__main__":
#    start_time = time.time()
#    pdb_ids_control = query_by_term_control("ribosome")
#    end_time = time.time()
#    print("Runtime for control: {:.2f} seconds".format(end_time - start_time))
   start_time = time.time()
   pdb_ids_exp = query_by_term("ribosome")  
   end_time = time.time()
   print("Runtime for experimental: {:.2f} seconds".format(end_time - start_time))
   print("this is output from pypdb:" + str(len(pdb_ids_exp)))
   with open('rcsb_pdb_ids_20231127152439.txt', 'r') as f:
       pdb_ids = f.read().split(',')
   with open('nakb_prna_ids.txt', 'r') as f:
       nakb_ids = f.read().split(',')
   print("this is output from rcsb:" + str(len(pdb_ids)))
   print("this is matches bw rcsb and our nakb dataset: " + str(len([entry.upper() for entry in nakb_ids if entry.upper() in set(pdb_ids)])))




