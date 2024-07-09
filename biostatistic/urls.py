from django.urls import path
from django.views.generic import TemplateView
from .views import create_table

app_name = "biostatistic"

urlpatterns = [
    path('anasayfa/', TemplateView.as_view(template_name="biostatistic/home.html"), name="home"),
    path('table-create/', create_table, name="table_create"),
]
