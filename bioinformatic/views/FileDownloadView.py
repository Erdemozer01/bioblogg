import os, mimetypes
from django.http.response import HttpResponseRedirect, HttpResponse
from pathlib import Path
from bioinformatic.models.bioinformatic import BioinformaticModel

BASE_DIR = Path(__file__).resolve().parent.parent.parent


def download_file(request):
    obj = BioinformaticModel.objects.filter(user=request.user)
    format = [j.writing_file_format for j in obj][0]
    filename = f'{request.user}_{format}_file.{format}'
    filepath = os.path.join(BASE_DIR, 'media', 'laboratory', f'{request.user}_{format}_file.{format}')
    path = open(filepath, 'r')
    mime_type, _ = mimetypes.guess_type(filepath)
    response = HttpResponse(path, content_type=mime_type)
    response['Content-Disposition'] = "attachment; filename=%s" % filename
    try:
        return response
    finally:
        path.close()
        os.remove(filepath)
        obj.delete()
