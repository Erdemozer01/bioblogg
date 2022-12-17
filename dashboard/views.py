from django.views import generic
from accounts.models import Profile
from django.contrib import messages
from django.shortcuts import redirect
from django.conf import settings
from blog.models import Comments


class DashboardView(generic.ListView):
    template_name = "dashboard/pages/dashboard.html"
    model = Profile

    def get(self, request, *args, **kwargs):
        if not request.user.is_staff:
            messages.error(request, "Yetkili girişi yapınız!")
            return redirect('%s?next=/blog/' % (settings.LOGIN_URL))
        return super().get(request, *args, **kwargs)

    def get_context_data(self, *, object_list=None, **kwargs):
        context = super().get_context_data(**kwargs)
        context['comments'] = Comments.objects.filter(post__author=self.request.user).order_by('-id')[:10]
        return context
