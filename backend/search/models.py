from django.db import models

class pypdbObject(models.Model):
    id = models.CharField(max_length=4, primary_key=True)
    description = models.CharField(max_length=1000)

    def __str__(self) -> str:
        return self.id + " " + self.description