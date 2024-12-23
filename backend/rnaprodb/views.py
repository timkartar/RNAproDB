from django.shortcuts import render
from rest_framework import viewsets
from .serializers import RNASerializer
from .models import RNA
from django.http import JsonResponse
import subprocess
import json
import os
from pypdb import get_info
# Create your views here.

class RNAView(viewsets.ModelViewSet):
    serializer_class = RNASerializer
    queryset = RNA.objects.all()

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

        if subgraph_nodes:
            script_path = "./get_subgraph.py"
            result = subprocess.run(["python", script_path, pdbid, subgraph_nodes, algorithm], capture_output=True, text=True)

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
                return JsonResponse({"message": "Error: No valid JSON found in the script's output."})
            
            try:
                json_output = json.loads(json_output)
                json_output['tooLarge'] = False
            except json.JSONDecodeError:
                return JsonResponse({"message": "Error decoding JSON output from script.", "error": errors})
            
            if result.returncode != 0:
                return JsonResponse({"message": "Error running script.", "error": errors})
        else: # full graph!
            RUN_RNAVIS_FLAG = True 
            if(RUN_RNAVIS_FLAG and algorithm == "pca"):
                script_path = "./rna_vis.py"
                result = subprocess.run(["python", script_path, pdbid], capture_output=True, text=True)
                
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
                    return JsonResponse({"message": "Error: No valid JSON found in the script's output."})
                
                try:
                    json_output = json.loads(json_output)
                except json.JSONDecodeError:
                    return JsonResponse({"message": "Error decoding JSON output from script.", "error": errors})
                
                if result.returncode != 0:
                    return JsonResponse({"message": "Error running script.", "error": errors})
               
                with open("./output/{}_{}_graph.json".format(pdbid, algorithm), 'r') as json_file:
                    json_output = json.load(json_file)  
            
            # DO NOT RUN RNA_VIS, FILES ALREADY THERE, meant for production
            else:
                with open("./output/{}_{}_graph.json".format(pdbid, algorithm), 'r') as json_file:
                    json_output = json.load(json_file)       
        # Use PyPDB to get title
        pdb_info = get_info(pdbid)

        TOO_LARGE = False
        if(json_output['tooLarge']):
            TOO_LARGE = True
        
        response_data = None
        if(TOO_LARGE and not subgraph_nodes):
            # new_json_output = {"chainsList": json_output['chainsList'], "ss": json_output["ss"]}
            response_data = {
            'file_url': '/{}.tmp.cif.html'.format(pdbid), 
            "message": "Script ran successfully!",
            "title": pdb_info['citation'][0]['title'],
            'protein_name': (pdb_info['struct']['title']),#.capitalize().replace('rna', 'RNA'),
            'tooLarge': True,
            "output": json_output,
            # "pdb_info": pdb_info
            }
        else:
            response_data = {
                'file_url': '/{}.tmp.cif.html'.format(pdbid), 
                "message": "Script ran successfully!",
                "title": pdb_info['citation'][0]['title'],
                'protein_name': (pdb_info['struct']['title']),#.capitalize().replace('rna', 'RNA'),
                'tooLarge': False,
                "output": json_output,  # Use the parsed JSON data here
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
    with open("./output/{}_{}_graph.json".format(pdbid, algorithm), 'r') as json_file:
        json_output = json.load(json_file)
    return JsonResponse(json_output)