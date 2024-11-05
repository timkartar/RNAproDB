from django.shortcuts import render
from rest_framework import viewsets
from .serializers import RNASerializer
from .models import RNA
from django.http import JsonResponse
import subprocess
import json
import os
from pypdb import get_info
import uuid
from django.views.decorators.csrf import ensure_csrf_cookie, csrf_exempt
from .make_table import makeTable
import re
# Create your views here.

# main_cwd = '/srv/www/rnaprodb/'
main_cwd = '/home/aricohen/Desktop/django-react-rnaprodb/'

temp_cwd = main_cwd + 'rnaprodb_dev'
MAX_FILE_SIZE = 50 * 1024 * 1024 # 50 MB
RUN_RNAVIS_FLAG = False

class RNAView(viewsets.ModelViewSet):
    serializer_class = RNASerializer
    queryset = RNA.objects.all()


def get_struct_list(request):
    struct_string = ""
    file_path = os.path.join(temp_cwd, "nakb_prna_ids.txt")
    
    try:
        with open(file_path, "r") as ids_file:
            struct_string = ids_file.read()  # Read the entire file content into struct_string
    except FileNotFoundError:
        return JsonResponse({"error": "File not found."}, status=404)
    except Exception as e:
        return JsonResponse({"error": f"An error occurred: {str(e)}"}, status=500)

    return JsonResponse({"structures": struct_string})



# helper function to run the main script, called by uploading and just pulling a structure
# pdb is uuid for an uploaded structure
def run_rna_vis(algorithm, pdbid, isUpload=False):
    json_output = None
    if((RUN_RNAVIS_FLAG and algorithm == "pca") or isUpload):
        script_path = "./rna_vis.py"
        # result = subprocess.run(["/home/aricohen/anaconda3/envs/RNAproDB/bin/python", script_path, pdbid], capture_output=True, text=True, cwd=temp_cwd)
        result = subprocess.run([f"{temp_cwd}/run_rna_vis_server.sh", pdbid], capture_output=True, text=True, cwd=temp_cwd)

        # # You can capture the stdout or stderr for further use if needed
        output = result.stdout
        errors = result.stderr

        # # Split the output by line breaks
        # lines = output.strip().split('\n')

        # # Find the JSON line (starting from the end)
        # for line in lines:
        #     if line.startswith("{"):
        #         break
        # json_output = line

        # if not json_output:
        #     return {"message": "Error: No valid JSON found in the script's output.", "output": output, "error": errors}
        
        try:
            with open("{}/output/{}_{}_graph.json".format(temp_cwd, pdbid, algorithm), 'r') as json_file:
                json_output = json.load(json_file)  
        except Exception as e:
            return {"message": "Error decoding JSON output from script.", "error": str(e), "output": output}
        
        # if result.returncode != 0:
        #     return {"message": "Error running script.", "error": errors}
        with open("{}/output/{}_{}_graph.json".format(temp_cwd, pdbid, algorithm), 'r') as json_file:
            json_output = json.load(json_file)  
            return json_output
    # DO NOT RUN RNA_VIS, FILES ALREADY THERE, meant for production
    else:
        with open("{}/output/{}_{}_graph.json".format(temp_cwd, pdbid, algorithm), 'r') as json_file:
            json_output = json.load(json_file)
            return json_output
"""
def run_electrostatics(request):
    if request.method == "GET":
        pdbid = request.GET.get('pdbid')
        
        # Check if pdbid is provided
        if not pdbid:
            return JsonResponse({"message": "Missing pdbid parameter."}, status=400)
        
        # Sanitize pdbid to allow only alphanumeric characters, dashes, and periods
        if not re.match(r'^[\w\-.]+$', pdbid):
            return JsonResponse({"message": "Invalid pdbid parameter."}, status=400)

        # Define the script path
        electro_path = os.path.join(temp_cwd, "electrostatics")
        script_path = os.path.join(electro_path, "process_upload.sh")

        # Check if the script exists
        if not os.path.isfile(script_path):
            return JsonResponse({"message": "Script not found."}, status=500)

        try:
            # Run the script synchronously
            result = subprocess.run(
                [script_path, pdbid],
                cwd=electro_path,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True  # Ensure stdout and stderr are captured as strings
            )
            
            # Check if there was an error in execution
            if result.returncode != 0:
                # Log and return the stderr output if there was an error
                print("Script stderr:", result.stderr)
                return JsonResponse({"message": f"Error running script: {result.stderr}"}, status=500)

            # Log the stdout output for successful execution
            print("Script stdout:", result.stdout)
            return JsonResponse({"message": f"Electrostatics process completed for {pdbid}.", "output": result.stdout}, status=200)

        except Exception as e:
            # Catch any exception and return an error response
            return JsonResponse({"message": f"Error starting process: {str(e)}"}, status=500)
"""
def run_electrostatics(request):
    if request.method == "GET":
        pdbid = request.GET.get('pdbid')
        
        # Check if pdbid is provided
        if not pdbid:
            return JsonResponse({"message": "Missing pdbid parameter."}, status=400)
        
        # Sanitize pdbid to allow only alphanumeric characters, dashes, and periods
        if not re.match(r'^[\w\-.]+$', pdbid):
            return JsonResponse({"message": "Invalid pdbid parameter."}, status=400)

        # Define the script path
        electro_path = os.path.join(temp_cwd, "electrostatics")
        script_path = os.path.join(electro_path, "process_upload.sh")

        # Check if the script exists
        if not os.path.isfile(script_path):
            return JsonResponse({"message": "Script not found."}, status=500)

        try:
            # Run the script asynchronously
            process = subprocess.Popen(
                [script_path, pdbid],
                cwd=electro_path,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )
            
            # Optionally, you could log the PID of the background process
            print(f"Started electrostatics process for {pdbid} with PID: {process.pid}")
            
            # Return response indicating the process has started
            return JsonResponse({"message": f"Electrostatics process started for {pdbid}."}, status=202)

        except Exception as e:
            # Catch any exception and return an error response
            return JsonResponse({"message": f"Error starting process: {str(e)}"}, status=500)

def run_script(request):
    # Ensure it's a GET request (although this will be the case by default for this route)
    if request.method == "GET":
        # Define the path to your script
        pdbid = request.GET.get('pdbid')
        subgraph_nodes = request.GET.get('subgraph')
        algorithm = request.GET.get('algorithm')
        isFirst = request.GET.get('isFirst', 'false').lower() == 'true'  # Correctly handle the 'isFirst' flag as boolean

        if not pdbid:
            return JsonResponse({"message": "Missing pdbid parameter."}, status=400)

        result = None
        json_output = None
        output_dir = "./output"
        table = None
        if subgraph_nodes:
            # script_path = "./get_subgraph.py"
            # result = subprocess.run(["/home/aricohen/anaconda3/envs/RNAproDB/bin/python", script_path, pdbid, subgraph_nodes, algorithm], capture_output=True, text=True, cwd=temp_cwd)
            result = subprocess.run([f"{temp_cwd}/run_subgraph_server.sh", pdbid, subgraph_nodes, algorithm], capture_output=True, text=True, cwd=temp_cwd)
            # result = subprocess.run([f"{temp_cwd}/run_subgraph_local.sh", pdbid, subgraph_nodes, algorithm], capture_output=True, text=True, cwd=temp_cwd)

            # You can capture the stdout or stderr for further use if needed
            output = result.stdout
            errors = result.stderr

            # Split the output by line breaks
            lines = output.strip().split('\n')

            # Find the JSON line (starting from the end)
            for line in lines:
                if line.startswith("'\"{"):
                    break
            json_output = line

            if not json_output:
                return JsonResponse({"message": "Error: No valid JSON found in the script's output.", "output": output, "errors": errors})
            
            try:
                json_output = json.loads(json_output)
                json_output['tooLarge'] = False
            except json.JSONDecodeError:
                return JsonResponse({"message": "Error decoding JSON output from script.", "error": errors})
            
            try:
                table = makeTable(json_output)
            except Exception as e:
                return JsonResponse({'error': str(e)}, status=400)

            if result.returncode != 0:
                return JsonResponse({"message": "Error running script.", "error": errors})
        else: # full graph!
            json_output = run_rna_vis(algorithm, pdbid)
            try:
                table = makeTable(json_output)
            except Exception as e:
                return JsonResponse({'error': str(e)}, status=400)
            if 'error' in json_output:
                return JsonResponse(json_output, status=400)
        # Use PyPDB to get title
        pdb_info = get_info(pdbid) 
        temp_title = "Uploaded Structure"
        temp_protein_name = ("N/A")

        # IF NOT UPLOADED STRUCTURE (longer PDB ID), get info
        if len(pdbid) <= 6:
            temp_title = pdb_info['citation'][0]['title']
            temp_protein_name = (pdb_info['struct']['title'])
        TOO_LARGE = False
        if(json_output['tooLarge']):
            TOO_LARGE = True
        
        response_data = None

        if not subgraph_nodes:
            # new_json_output = {"chainsList": json_output['chainsList'], "ss": json_output["ss"]}
            response_data = {
            'file_url': '/{}.tmp.cif.html'.format(pdbid), 
            "message": "Script ran successfully!",
            "title": temp_title,
            'protein_name': temp_protein_name,#.capitalize().replace('rna', 'RNA'),
            'tooLarge': False,
            "table": table,
            "output": json_output,
            # "pdb_info": pdb_info
            }
        else:
            response_data = {
                'file_url': '/{}.tmp.cif.html'.format(pdbid), 
                "message": "Script ran successfully!",
                "title": temp_title,
                'protein_name': temp_protein_name,#.capitalize().replace('rna', 'RNA'),
                'tooLarge': False,
                "output": json_output,  # Use the parsed JSON data here
                "table": table
                # "pdb_info": pdb_info
            }
        return JsonResponse(response_data)
        
    return JsonResponse({"message": "Invalid request method."}, status=405)

def get_struct_info(request):
    pdbid = request.GET.get('pdbid')
    pdb_info = get_info(pdbid)
    return JsonResponse(pdb_info)

def download_json(request):
    pdbid = request.GET.get('pdbid')
    algorithm = request.GET.get('algorithm')
    with open("{}/output/{}_{}_graph.json".format(temp_cwd, pdbid, algorithm), 'r') as json_file:
        json_output = json.load(json_file)
    return JsonResponse(json_output)

# for now no async
@csrf_exempt
def handle_upload(request):
    file = request.FILES.get('file')

    if not file:
        return JsonResponse({"error": "No file provided"}, status=400)
    
    if file.size > MAX_FILE_SIZE:
        return JsonResponse({'error': 'File size exceeds the allowed limit'}, status=400)

    unique_id = "upload-" + str(uuid.uuid4())
    file_extension = os.path.splitext(file.name)[1]
    
    if file_extension.lower() not in ['.cif']:
        return JsonResponse({'error': 'Invalid file type'}, status=400)

    # Define the path where the file will be saved
    file_path = os.path.join(f'{main_cwd}/rnaprodb_dev/output/cifs', f'{unique_id}-assembly1{file_extension}')

    # Write the file to the disk
    with open(file_path, 'wb+') as destination:
        for chunk in file.chunks():
            destination.write(chunk)
    
    json_output = run_rna_vis('pca', unique_id, isUpload=True)
    if json_output:
        if 'error' in json_output: # if dictionary did not work!
            return JsonResponse({"message": "error processing your file", "error": "error processing your file", "output": str(json_output)}, status=400)
        response_data = {
                "message": "Script ran successfully!",
                "id": unique_id,
        }
        return JsonResponse(response_data)
    else:
        return JsonResponse({'error': 'Upload failed script no run.', 'json_output': json_output}, status=400)
