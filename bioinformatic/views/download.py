import mimetypes
import os
from django.http.response import HttpResponse
from pathlib import Path
from bioinformatic.models.bioinformatic import BioinformaticModel

BASE_DIR = Path(__file__).resolve().parent.parent.parent


def download_file(request):
    obj = BioinformaticModel.objects.filter(user=request.user)
    format = [j.writing_file_format for j in obj][0]
    # Define text file name
    filename = f'{request.user}_{format}_file.{format}'
    # Define the full file path
    filepath = os.path.join(BASE_DIR, 'media', f'{request.user}_{format}_file.{format}')
    # Open the file for reading content
    path = open(filepath, 'r')
    # Set the mime type
    mime_type, _ = mimetypes.guess_type(filepath)
    # Set the return value of the HttpResponse
    response = HttpResponse(path, content_type=mime_type)
    # Set the HTTP header for sending to browser
    response['Content-Disposition'] = "attachment; filename=%s" % filename
    try:
        # Return the response value
        return response
    finally:
        path.close()
        os.remove(filepath)
        obj.delete()

