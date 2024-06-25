import pypdb
from pypdb.clients.search.search_client import perform_search, perform_search_with_graph, ReturnType, QueryGroup, LogicalOperator
from pypdb.clients.search.operators import text_operators
import sys


def query_by_term(term: str, queries: list) -> list:
   search_query = text_operators.DefaultOperator(value = term)
   queries.append(search_query)

def query_by_resolution(min_res:float, max_res:float, queries: list):
   resolution_query = text_operators.RangeOperator(
         attribute="rcsb_entry_info.resolution_combined",
         from_value=min_res,
         to_value=max_res,
         include_lower=True,
         include_upper=True,
         negation=False)
   queries.append(resolution_query)

def query_by_NA(min_NA:float, max_NA:float, queries: list):
   NA_query = text_operators.RangeOperator(
         attribute="rcsb_assembly_info.polymer_entity_count_RNA",
         from_value=min_NA,
         to_value=max_NA,
         include_lower=True,
         include_upper=True,
         negation=False)
   queries.append(NA_query)

def query_by_protein(min_protein:float, max_protein:float, queries: list):
   protein_query = text_operators.RangeOperator(
      attribute="rcsb_assembly_info.polymer_entity_count_protein",
      from_value=min_protein,
      to_value=max_protein,
      include_lower=True,
      include_upper=True,
      negation=False)
   queries.append(protein_query)

def query_by_experimental_modality(modality: list, queries: list):
   modality_query = text_operators.InOperator(
      attribute="rcsb_entry_info.experimental_method",
      values=modality)
   queries.append(modality_query)

def query_by_year(min_year:int, max_year:int, queries: list):
   year_query = text_operators.RangeOperator(
      attribute="rcsb_primary_citation.year",
      from_value=min_year,
      to_value=max_year,
      include_lower=True,
      include_upper=True,
      negation=False)
   queries.append(year_query)

def search(additional_queries: list) -> list:
   queries = []
   queries.append(pypdb.text_operators.ComparisonOperator(value = 0, attribute = "rcsb_assembly_info.polymer_entity_count_RNA", comparison_type = pypdb.text_operators.ComparisonType.GREATER))
   queries.append(pypdb.text_operators.ComparisonOperator(value = 0, attribute = "rcsb_assembly_info.polymer_entity_count_protein", comparison_type = pypdb.text_operators.ComparisonType.GREATER))
   for i in additional_queries:
      if i[0] == "term":
         query_by_term(i[1], queries)
      elif i[0] == "resolution":
         query_by_resolution(i[1][0], i[1][1], queries)
      elif i[0] == "NA":
         query_by_NA(i[1][0], i[1][1], queries)
      elif i[0] == "protein":
         query_by_protein(i[1][0], i[1][1], queries)
      elif i[0] == "experimental_modality":
         query_by_experimental_modality(i[1], queries)
      elif i[0] == "year":
         query_by_year(i[1][0], i[1][1], queries)

   results = perform_search_with_graph(
   query_object=QueryGroup(
      logical_operator=LogicalOperator.AND,
      queries=queries
   ),
   return_type=ReturnType.ENTRY
   )
   with open('nakb_prna_ids.txt', 'r') as f:
      nakb_ids = f.read().split(',')
   output_list = [entry.lower() for entry in nakb_ids if entry.upper() in set(results)]
   print(additional_queries)
   print("Output length: " + str(len(output_list)))
   return output_list

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
   pdb = sys.argv[1]
   print(pypdb.get_info(pdb))




