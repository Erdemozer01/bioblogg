from django.views import generic
from django.urls import reverse_lazy
from django.shortcuts import render, reverse, redirect
from django.views.generic import ListView
from blog.models import Posts
from django.contrib.auth.models import User
from django.contrib import messages
from .forms import UserRegistrationForm, UserEditForm, UserCreationForm, ProfileEditForm, UserProfileEditForm
from django.contrib.auth.mixins import LoginRequiredMixin
from .models import Profile
from django.contrib.auth.decorators import login_required
from django.contrib.auth.views import PasswordChangeView
from django.contrib.auth.forms import PasswordChangeForm
from django.views.generic import DetailView, CreateView


class UserRegister(generic.CreateView):
    template_name = "registration/sign-up.html"
    form_class = UserRegistrationForm
    success_url = reverse_lazy("accounts:register")

    def form_valid(self, form):
        messages.success(self.request, 'Başarılı Bir Şekilde Kayıt Oldunuz.')
        form.save()
        return super(UserRegister, self).form_valid(form)


