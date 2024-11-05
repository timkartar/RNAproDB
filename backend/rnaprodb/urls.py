from django.urls import path
from . import views

urlpatterns = [
    path('run-script/', views.run_script, name='run_script'),
    path('get_struct_info/', views.get_struct_info, name='get_struct_info'),
    path('download_json/', views.download_json, name='download_json'),
    path('handle_upload', views.handle_upload, name='handle_upload'),
    path('run-electrostatics', views.run_electrostatics, name='run_electrostatics'),
    path('get-struct-list', views.get_struct_list, name='get_struct_list'),
]