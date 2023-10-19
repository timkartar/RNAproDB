from django.shortcuts import render
from rest_framework import viewsets
from .serializers import RNASerializer
from .models import RNA
from django.http import JsonResponse
import subprocess
# Create your views here.

class RNAView(viewsets.ModelViewSet):
    serializer_class = RNASerializer
    queryset = RNA.objects.all()

def run_script(request):
    # Ensure it's a GET request (although this will be the case by default for this route)
    if request.method == "GET":
        # Define the path to your script
        pdbid = request.GET.get('pdbid')
        
        if not pdbid:
            return JsonResponse({"message": "Missing pdbid parameter."}, status=400)

        
        script_path = "/home/aricohen/Desktop/django-react-rnaprodb/rnaprodb_dev/rna_vis.py"
        result = subprocess.run(["python", script_path, pdbid], capture_output=True, text=True)

        # You can capture the stdout or stderr for further use if needed
        output = result.stdout
        errors = result.stderr

        # Return the result or handle errors as desired
        if result.returncode == 0:
            return JsonResponse({'file_url': '/{}.tmp.cif.html'.format(pdbid), "message": "Script ran successfully!", "output": output})
        else:
            return JsonResponse({"message": "Error running script.", "error": errors})
    
    return JsonResponse({"message": "Invalid request method."}, status=405)
