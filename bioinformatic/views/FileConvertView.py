from django.shortcuts import *
from pathlib import Path
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from django.contrib import messages
from bioinformatic.forms.file_convert import FileConvertForm
from bioinformatic.models.bioinformatic import BioinformaticModel
from django.utils.translation import gettext as _
from django.utils.translation import ngettext_lazy
BASE_DIR = Path(__file__).resolve().parent.parent.parent


def FileConvert(request):
    global in_file, obj, out_file
    if request.user.is_anonymous:
        from django.conf import settings
        messages.error(request, "Lütfen Giriş Yapınız")
        return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))

    form = FileConvertForm(request.POST or None, request.FILES or None)

    if BioinformaticModel.objects.filter(user=request.user, tool="DOSYA DÖNÜŞTÜRME").exists():
        BioinformaticModel.objects.filter(user=request.user, tool="DOSYA DÖNÜŞTÜRME").delete()

    if request.method == "POST":
        if form.is_valid():
            try:
                obj = BioinformaticModel.objects.create(
                    user=request.user,
                    tool="DOSYA DÖNÜŞTÜRME",
                    molecule=form.cleaned_data['molecule'],
                    reading_file_format=form.cleaned_data['reading_file_format'],
                    writing_file_format=form.cleaned_data['writing_file_format'],
                )

                in_file = obj.records_files.create(
                    file=request.FILES['file']
                )

                out_file = obj.records_files.create(
                    file=f"{request.user}_{obj.writing_file_format}_file.{obj.writing_file_format}"
                )

                SeqIO.convert(in_file.file.path, obj.reading_file_format, out_file.file.path, obj.writing_file_format,
                              obj.molecule)

            except ValueError as err:
                BioinformaticModel.objects.filter(user=request.user, tool="DOSYA DÖNÜŞTÜRME").delete()
                messages.error(request, str(err))
                return HttpResponseRedirect(request.get_full_path())

        return HttpResponseRedirect(reverse("bioinformatic:download"))

    return render(request, "bioinformatic/form.html", {'form': form, 'title': _("DOSYA DÖNÜŞTÜRME")})
