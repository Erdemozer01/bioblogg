from django.shortcuts import *
from pathlib import Path
from django.views import generic
from django.contrib import messages

BASE_DIR = Path(__file__).resolve().parent.parent.parent


def file_writing(request, user):
    if request.user.is_anonymous:
        from django.conf import settings
        messages.error(request, "Lütfen Giriş Yapınız")
        return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))

    return render(request, 'bioinformatic/form.html')
