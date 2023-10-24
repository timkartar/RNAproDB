import pickle
import csv

def pickle_to_csv(pickle_file_path: str, csv_name):

    with open(pickle_file_path, 'rb') as FH:
        data = pickle.load(FH)

    with open(csv_name, 'w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(['PDB ID', 'Count'])
        for pdb_id, count in data.items():
            csv_writer.writerow([pdb_id, count])

if __name__ == "__main__":
    #pickle_to_csv('hirad/rna_count_dict.pickle', 'RNA_nuc_counts.csv')
    pickle_to_csv('aa_counts_dict.pickle', 'aa_counts.csv')
