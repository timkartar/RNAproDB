from rest_framework.response import Response
from rest_framework.decorators import api_view
from rest_framework import status
from django.http import JsonResponse
from .models import pypdbObject
from .serializers import pypdbSerializer
from .pypdbSearch import *
import json
from requests.exceptions import JSONDecodeError

@api_view(['GET', 'POST'])
def search_view(request):
    if request.method == 'POST':
        try:
            data = json.loads(request.body)
            search_term = data.get('searchTerm')
            min_resolution = data.get('minResolution')
            max_resolution = data.get('maxResolution')
            min_NA = data.get('minNA')  
            max_NA = data.get('maxNA')
            min_protein = data.get('minProtein')
            max_protein = data.get('maxProtein')
            experimental_modality = data.get('experimentalModality')
            # Extract other parameters as needed

            additional_queries = []
            if search_term:
                additional_queries.append(['term', search_term])
            if min_resolution and max_resolution:
                additional_queries.append(['resolution', [min_resolution, max_resolution]])
            if min_NA and max_NA:
                additional_queries.append(['NA', [min_NA, max_NA]])
            if min_protein and max_protein:
                additional_queries.append(['protein', [min_protein, max_protein]])
            if experimental_modality:
               additional_queries.append(['experimental_modality', experimental_modality])
            desired_ids = search(additional_queries)
            pdbs = pypdbObject.objects.filter(id__in=desired_ids)
            serializer = pypdbSerializer(pdbs, many=True)
            return JsonResponse(serializer.data, safe=False)
        except JSONDecodeError as e:
            return Response({'error': 'Bad request', 'message': str(e)}, status=status.HTTP_400_BAD_REQUEST)
        except Exception as e:  # Catching other exceptions
            # Log the exception for debugging
            # log.error(f"Exception in search_view: {str(e)}")
            return Response({'error': 'Internal Server Error', 'message': str(e)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)

    elif request.method == 'GET':
        # Handle GET request if necessary
        # For example, return a default set of data or an error message
        return Response({"message": "GET requests are not supported on this endpoint"}, status=status.HTTP_405_METHOD_NOT_ALLOWED)

@api_view(['GET', 'POST'])
def pdb_list(request, format=None):
    if request.method == 'GET':
        desired_ids = query_by_term("ribosome")
        pdbs = pypdbObject.objects.filter(id__in=desired_ids)
        serializer = pypdbSerializer(pdbs, many=True)
        return Response(serializer.data)
    if request.method == 'POST':
        serializer = pypdbSerializer(data=request.data)
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data, status=status.HTTP_201_CREATED)

@api_view(['GET', 'PUT', 'DELETE'])
def pdb_detail(request, id, format=None):
    try:    
        pdb = pypdbObject.objects.get(pk=id)
    except pypdbObject.DoesNotExist:
        return Response(status=status.HTTP_404_NOT_FOUND)

    if request.method == 'GET':
        serializer = pypdbSerializer(pdb)
        return Response(serializer.data)
    elif request.method == 'PUT':
        serializer = pypdbSerializer(pdb, data=request.data)
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)
    elif request.method == 'DELETE':
        pdb.delete()
        return Response(status=status.HTTP_204_NO_CONTENT)