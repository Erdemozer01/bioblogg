from django.urls import path
from .views import (
    DashboardView
)

app_name = "dashboard"

urlpatterns = [
    path('<str:user>/', DashboardView.as_view(), name="anasayfa"),
]
