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

def query_by_NA_type(NA_types: list, queries: list):
    for NA_type in NA_types:
        if NA_type == "RNA (only)":
            NA_type_query = text_operators.ComparisonOperator(
                value=0, 
                attribute="rcsb_assembly_info.polymer_entity_count_RNA", 
                comparison_type=pypdb.text_operators.ComparisonType.GREATER
            )
            queries.append(NA_type_query)
        elif NA_type == "DNA (only)":
            NA_type_query = text_operators.ComparisonOperator(
                value=0, 
                attribute="rcsb_assembly_info.polymer_entity_count_DNA", 
                comparison_type=pypdb.text_operators.ComparisonType.GREATER
            )
            queries.append(NA_type_query)
        elif NA_type == "Hybrid":
            NA_type_query_RNA = text_operators.ComparisonOperator(
                value=0, 
                attribute="rcsb_assembly_info.polymer_entity_count_RNA", 
                comparison_type=pypdb.text_operators.ComparisonType.GREATER
            )
            NA_type_query_DNA = text_operators.ComparisonOperator(
                value=0, 
                attribute="rcsb_assembly_info.polymer_entity_count_DNA", 
                comparison_type=pypdb.text_operators.ComparisonType.GREATER
            )
            queries.append(NA_type_query_RNA)
            queries.append(NA_type_query_DNA)

def query_by_molecular_weight(min_mw:float, max_mw:float, queries: list): 
      mw_query = text_operators.RangeOperator(
         attribute="rcsb_entry_info.molecular_weight",
         from_value=min_mw,
         to_value=max_mw,
         include_lower=True,
         include_upper=True,
         negation=False)
      queries.append(mw_query)

def search(additional_queries: list) -> list:
    queries = []
    
    if additional_queries:
        for query_type, query_value in additional_queries:
            if query_type == "term":
                query_by_term(query_value, queries)
            elif query_type == "resolution":
                query_by_resolution(query_value[0], query_value[1], queries)
            elif query_type == "NA":
                query_by_NA(query_value[0], query_value[1], queries)
            elif query_type == "protein":
                query_by_protein(query_value[0], query_value[1], queries)
            elif query_type == "experimental_modality":
                query_by_experimental_modality(query_value, queries)
            elif query_type == "year":
                query_by_year(query_value[0], query_value[1], queries)
            elif query_type == "NA_type":
                query_by_NA_type(query_value, queries)
            elif query_type == "molecular_weight":
                query_by_molecular_weight(query_value[0], query_value[1], queries)

        results = perform_search_with_graph(
            query_object=QueryGroup(
                logical_operator=LogicalOperator.AND,
                queries=queries
            ),
            return_type=ReturnType.ENTRY
        )
        result_set = {entry.upper() for entry in results}
    else:
        result_set = set()

    # Read the file only once
    with open('nakb_prna_ids.txt', 'r') as f:
        nakb_ids = {entry.strip().upper() for entry in f.read().split(',')}
    
    if not additional_queries:
        output_set = nakb_ids
    else:
        output_set = nakb_ids.intersection(result_set)
    
    output_list = [entry.lower() for entry in output_set]
    
    print(additional_queries)
    print("Output length: " + str(len(output_list)))
    return output_list

# Example usage:
if __name__ == "__main__":
    import sys
    pdb = sys.argv[1]
    print(pypdb.get_info(pdb))






