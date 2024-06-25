from rest_framework import serializers
from .models import pypdbObject

class pypdbSerializer(serializers.ModelSerializer):
    class Meta:
        model = pypdbObject
        fields = ['id', 'authors', 'title', 'year', 'doi', 'pubmed']