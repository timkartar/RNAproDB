from rest_framework import serializers
from .models import RNA

class RNASerializer(serializers.ModelSerializer):
    class Meta:
        model = RNA
        fields = ('id', 'name', 'chain', 'position')