from rest_framework.response import Response
from rest_framework.decorators import api_view
from rest_framework import status
from django.http import JsonResponse
from .models import pypdbObject
from .serializers import pypdbSerializer
from .pypdbSearch import *

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