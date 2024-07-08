from django.urls import path
from django.views.generic import TemplateView

app_name = "biostatistic"

urlpatterns = [
    path('anasayfa/', TemplateView.as_view(template_name="biostatistic/home.html"), name="home"),
]
