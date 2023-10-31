from django.shortcuts import render
from rest_framework import viewsets
from .serializers import RNASerializer
from .models import RNA
from django.http import JsonResponse
import subprocess
import json
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

        if not pdbid:
            return JsonResponse({"message": "Missing pdbid parameter."}, status=400)

        # script_path = "/home/raktim/rnaprodb/rnaprodb/rna_vis.py"
        # script_path = "./rna_vis.py"
        result = None
        json_output = None
        if subgraph_nodes:
            script_path = "./get_subgraph.py"
            result = subprocess.run(["python", script_path, pdbid, subgraph_nodes], capture_output=True, text=True)

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
        else: # full graph!
            # script_path = "./rna_vis.py"
            with open("./output/{}_graph.json".format(pdbid), 'r') as json_file:
                json_output = json.load(json_file)

        # Use PyPDB to get title
        pdb_info = get_info(pdbid)

        response_data = {
            'file_url': '/{}.tmp.cif.html'.format(pdbid), 
            "message": "Script ran successfully!",
            "title": pdb_info['citation'][0]['title'],
            'protein_name': (pdb_info['struct']['title']),#.capitalize().replace('rna', 'RNA'),
            "output": json_output  # Use the parsed JSON data here
        }
        return JsonResponse(response_data)
        
    return JsonResponse({"message": "Invalid request method."}, status=405)
