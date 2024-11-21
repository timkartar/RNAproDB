from django.db import models

class pypdbObject(models.Model):
    id = models.CharField(max_length=4, primary_key=True)
    authors = models.CharField(max_length=1000, null=True)
    title = models.CharField(max_length=1000, null=True)
    year = models.IntegerField(null=True)
    doi = models.CharField(max_length=100, null=True)
    pubmed = models.IntegerField(null=True)
    is_rna_protein = models.BooleanField(null=False)

    def __str__(self) -> str:
        return self.id + ": " + self.title