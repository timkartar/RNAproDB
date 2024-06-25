from django.urls import path
from . import views
from rest_framework.urlpatterns import format_suffix_patterns

urlpatterns = [
    path('pypdb/', views.search_view),
    path('pypdb/<str:id>/', views.pdb_detail),
]

urlpatterns = format_suffix_patterns(urlpatterns)