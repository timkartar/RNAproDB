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
# Create your views here.

main_cwd = '/srv/www/rnaprodb/'
# main_cwd = '/home/aricohen/Desktop/django-react-rnaprodb/'

temp_cwd = main_cwd + 'rnaprodb_dev'
MAX_FILE_SIZE = 50 * 1024 * 1024 # 50 MB
RUN_RNAVIS_FLAG = False

class RNAView(viewsets.ModelViewSet):
    serializer_class = RNASerializer
    queryset = RNA.objects.all()

# helper function to run the main script, called by uploading and just pulling a structure
# pdb is uuid for an uploaded structure
def run_rna_vis(algorithm, pdbid, isUpload=False):
    json_output = None
    if((RUN_RNAVIS_FLAG and algorithm == "pca") or isUpload):
        # script_path = "./rna_vis.py"
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
            json_output = json.loads(f"{temp_cwd}/output/upload-{pdbid}_pca_graph.json")
        except json.JSONDecodeError:
            return {"message": "Error decoding JSON output from script.", "error": errors, "output": output}
        
        if result.returncode != 0:
            return {"message": "Error running script.", "error": errors}
        with open("{}/output/{}_{}_graph.json".format(temp_cwd, pdbid, algorithm), 'r') as json_file:
            json_output = json.load(json_file)  
            return json_output
    # DO NOT RUN RNA_VIS, FILES ALREADY THERE, meant for production
    else:
        with open("{}/output/{}_{}_graph.json".format(temp_cwd, pdbid, algorithm), 'r') as json_file:
            json_output = json.load(json_file)
            return json_output
# refactor to work with uploads
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
    file_path = os.path.join(f'{main_cwd}/rnaprodb_frontend/build/cifs', f'{unique_id}-assembly1{file_extension}')

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
