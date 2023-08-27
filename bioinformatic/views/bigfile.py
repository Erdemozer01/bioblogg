from django.shortcuts import *
from django.conf import settings

def big_file_upload(request):
    # Define Django project base directory
    BASE_DIR = settings.BIG_FILE_STORAGE
    # Define text file name
    filename = request.FILES['file']