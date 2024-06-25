from django.urls import path
from . import views

urlpatterns = [
    path('run-script/', views.run_script, name='run_script'),
    path('get_struct_info/', views.get_struct_info, name='get_struct_info'),
]