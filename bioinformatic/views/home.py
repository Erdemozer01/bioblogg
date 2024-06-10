from django.shortcuts import render


def bioinformatic_home(request):
    return render(request, "bioinformatic/home.html")
