from django.contrib.auth.models import User
from django.views import generic
from django.contrib import messages
from .forms import UserRegistrationForm
from django.urls import reverse_lazy
from django.core.mail import send_mail

class UserRegister(generic.CreateView):
    template_name = "registration/sign-up.html"
    form_class = UserRegistrationForm
    success_url = reverse_lazy("accounts:register")

    def form_valid(self, form):
        user = form.cleaned_data['username']
        user_count = User.objects.all().count()
        messages.success(self.request, 'Başarılı Bir Şekilde Kayıt Oldunuz.')
        send_mail(
            "Yeni üye",
            f"{user} kullanıcı adıyla siteye üye oldu. \n Sitenin toplam üye sayısı {user_count} ulaştı.",
            "bioblogdestek@gmail.com",
            ["ozer246@gmail.com"],
            fail_silently=False
        )
        form.save()
        return super(UserRegister, self).form_valid(form)
