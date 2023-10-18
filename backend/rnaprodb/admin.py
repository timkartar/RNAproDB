from django.contrib import admin
from .models import RNA

class RnaProDbAdmin(admin.ModelAdmin):
    list_display = ('name', 'chain', 'position')
# Register your models here.

admin.site.register(RNA,RnaProDbAdmin)