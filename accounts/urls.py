from django.urls import path
from accounts.views import UserRegister

app_name = "accounts"

urlpatterns = [
    path('register/', UserRegister.as_view(), name="register"),

]
