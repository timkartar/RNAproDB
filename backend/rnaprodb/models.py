from django.db import models

# Create your models here.
class RNA(models.Model):
    name = models.CharField(max_length=120)
    chain = models.CharField(max_length=120)
    position = models.IntegerField()
    def _str_(self):
        return str(self.name) + ":" + str(self.chain) + ":" + str(self.position)