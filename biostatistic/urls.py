from django.urls import path
from django.views.generic import TemplateView
from .views import files_table

app_name = "biostatistic"

urlpatterns = [
    path('anasayfa/', TemplateView.as_view(template_name="biostatistic/home.html"), name="home"),
    path('files-table/', files_table, name="files_table"),
]
