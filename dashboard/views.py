from django.views import generic
from accounts.models import Profile
from django.contrib import messages
from django.shortcuts import redirect
from django.conf import settings


class DashboardView(generic.ListView):
    template_name = "dashboard/pages/dashboard.html"
    model = Profile

    def get(self, request, *args, **kwargs):
        if request.user.is_anonymous:
            messages.error(request, "Lütfen Giriş Yapınız")
            return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))
        return super().get(request, *args, **kwargs)
