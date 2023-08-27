import uuid
import mimetypes
from pathlib import Path
import os
from django.shortcuts import *
from bioinformatic.forms.convert import FileConvertForm, FileConvertModelForm
from Bio import SeqIO
from django.views import generic
from django.http import FileResponse

BASE_DIR = Path(__file__).resolve().parent.parent.parent


def handle_uploaded_file(f):
    with open(os.path.join(BASE_DIR, "bioinformatic", "files", f"{f}"), 'wb+') as destination:
        for chunk in f.chunks():
            destination.write(chunk)


def file_convert(request):
    form = FileConvertForm(request.POST or None, request.FILES or None)
    request.user.id = uuid.uuid4().int
    if request.method == "POST":
        if form.is_valid():
            handle_uploaded_file(request.FILES['in_file'])
            in_format = form.cleaned_data['in_format']
            out_format = form.cleaned_data['out_format']
            in_file_path = os.path.join(BASE_DIR, "bioinformatic", "files", f"{request.FILES['in_file']}")
            out_file_path = os.path.join(BASE_DIR, "bioinformatic", "files",
                                         "{}_converted.{}".format(in_format, out_format))
            molecule_type = form.cleaned_data['molecule_type']
            try:
                SeqIO.convert(
                    in_file=in_file_path,
                    in_format=in_format,
                    out_file=out_file_path,
                    out_format=out_format,
                    molecule_type=molecule_type
                )
                path = open(out_file_path, 'rb')
                # Set the mime type
                mime_type, _ = mimetypes.guess_type(out_file_path)
                # Set the return value of the HttpResponse
                response = HttpResponse(path, content_type=mime_type)
                # Set the HTTP header for sending to browser
                response['Content-Disposition'] = "attachment; filename=%s" % "{}_converted_{}.{}".format(in_format,
                                                                                                          str(request.user.id)[
                                                                                                          :5],
                                                                                                          out_format)

                handle = open(in_file_path)
                handle.close()
                os.remove(in_file_path)

                return response

            finally:
                handle = open(out_file_path)
                handle.close()
                os.remove(out_file_path)

    return render(request, "bioinformatic/convert/input.html", {'form': form, 'bre': 'Dosya Dönüştürme'})


from django.contrib.auth.decorators import permission_required, login_required
from bioinformatic.models import FileFormat
from django.contrib.auth import settings


@login_required
@permission_required('bioinformatic.add_file_format', login_url=settings.LOGIN_URL, raise_exception=True)
def add_file_format(request):
    form = FileConvertModelForm(request.POST or None)
    if request.method == "POST":
        if form.is_valid():
            file_format = form.cleaned_data['name']
            if FileFormat.objects.filter(name=file_format).exists():
                pass
            else:
                form.save()
            return redirect('bioinformatic:file_convert')
        else:
            form = FileConvertModelForm()
    return render(request, "bioinformatic/convert/add_format.html", {'form': form, 'bre': 'Dosya Formatı Ekle'})


@login_required
@permission_required('bioinformatic.add_file_format', login_url=settings.LOGIN_URL, raise_exception=True)
def file_format_delete(request, pk):
    obj = FileFormat.objects.get(pk=pk)
    obj.delete()
    if FileFormat.DoesNotExist:
        return redirect('bioinformatic:file_convert')
    return redirect('bioinformatic:edit_file_format')


@login_required
@permission_required('bioinformatic.add_file_format', login_url=settings.LOGIN_URL, raise_exception=True)
def file_formats(request):
    object_list = FileFormat.objects.all()
    return render(request, "bioinformatic/convert/add_format.html",
                  {'object_list': object_list, 'bre': 'Dosya Formatı Düzenle'})
