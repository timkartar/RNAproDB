# get list of relevant IDS
# cross reference with existing ones/already considered!
    #include the date for these
# download the first biological assembly
# run rna_vis
# update the list

# IGNORE SEARCH FOR NOW

import subprocess
import os
import json
import download_ids

home =  os.path.dirname(os.path.abspath(__file__))
backend =  os.path.dirname(os.path.abspath(__file__))
frontend = backend + "/../rnaprodb_frontend/"

pdb_path = frontend + "build/cifs/"

"""
Returns true if worked, false if not
"""
def run_rna_vis(pdbid):    
    # Construct the command as a list. This runs `python rna_vis.py pdbid`
    command = ["python", "rna_vis.py", pdbid]

    try:
        # Run the command and capture the output
        result = subprocess.run(command, capture_output=True, text=True, check=True)

    except subprocess.CalledProcessError as e:
        # Handle errors in the subprocess
        print(f"An error occurred: {e}")
        print(f"Return code: {e.returncode}")
        print(f"Output: {e.output}")
        print(f"Error: {e.stderr}")

    # Construct the file path
    file_path = '{}/output/{}_pca_graph.json'.format(home, pdbid)

    # Check if the file exists
    if os.path.exists(file_path):
        try:
            # Try to open the file and parse its content as JSON
            with open(file_path, 'r') as outfile:
                data = json.load(outfile)
                print("Valid JSON file")
                return True
        except json.JSONDecodeError:
            print("Error: File contains invalid JSON.")
            return False
        except Exception as e:
            print(f"An error occurred: {e}")
            return False
    else:
        print(f"File does not exist: {file_path}")
        return False
    return False

if __name__ == "__main__":
    potential_ids = set(download_ids.main())
    seen_ids = set()
    with open('ids_considered.txt', 'r') as ids_file:
        for cur_id in ids_file:
            seen_ids.add(cur_id.strip())
    
    # New IDS not common to both
    new_ids = potential_ids.symmetric_difference(seen_ids)

    to_run_electrostatics = []

    with open('ids_considered.txt', 'a') as ids_file:
        for new_id in new_ids:
            ran_success = run_rna_vis(new_id)
            if not ran_success:
                print('Failed', new_id)
            else:
                print('Succeeded', new_id)
                print('Adding to electrostatics queue')
                to_run_electrostatics.append(new_id)
            ids_file.write("\n" + new_id)
 
    # Run search update

    # Run electrostatics
    electrostatics_dir = backend + "/electrostatics"
    for id in to_run_electrostatics:
        try:
            # Run the script with the specified ID in the correct directory
            subprocess.run(
                ["./process_upload.sh", id],
                cwd=electrostatics_dir,
                check=True,          # Raises an error if the command fails
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True             # Ensures stdout/stderr are strings
            )
            print(f"Electrostatics process completed for {id}")
        except subprocess.CalledProcessError as e:
            # Handle any errors that occur during the process
            print(f"Error running electrostatics for {id}: {e.stderr}")
