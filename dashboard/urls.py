from django.urls import path
from .views import (
    DashboardView,
    BlogDashBoardView
)

app_name = "dashboard"

urlpatterns = [
    path('<str:user>/', DashboardView.as_view(), name="anasayfa"),
    path('my-blog-dashboard/<pk>/<user>/', BlogDashBoardView.as_view(), name="my_blog_dashboard"),
]
