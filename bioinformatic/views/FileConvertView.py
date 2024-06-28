from django.shortcuts import *
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from django.contrib import messages
from bioinformatic.forms.file_convert import FileConvertForm
from bioinformatic.models.bioinformatic import BioinformaticModel
from django.utils.translation import gettext as _
from django.http import FileResponse

BASE_DIR = Path(__file__).resolve().parent.parent.parent


def FileConvert(request):
    global in_handle, obj, out_handle
    if request.user.is_anonymous:
        from django.conf import settings
        messages.error(request, "Lütfen Giriş Yapınız")
        return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))

    form = FileConvertForm(request.POST or None, request.FILES or None)

    if BioinformaticModel.objects.filter(user=request.user, tool="DOSYA DÖNÜŞTÜRME").exists():
        BioinformaticModel.objects.filter(user=request.user, tool="DOSYA DÖNÜŞTÜRME").delete()

    if request.method == "POST":
        if form.is_valid():

            obj = BioinformaticModel.objects.create(
                user=request.user,
                tool="DOSYA DÖNÜŞTÜRME",
                molecule=form.cleaned_data['molecule'],
                reading_file_format=form.cleaned_data['reading_file_format'],
                writing_file_format=form.cleaned_data['writing_file_format'],
            )

            in_handle = obj.records_files.create(
                file=request.FILES['file']
            )

            out_handle = obj.records_files.create(
                file=f"{request.user}_{obj.writing_file_format}_file.{obj.writing_file_format}"
            )

            records = SeqIO.parse(in_handle.file.path, obj.reading_file_format)

            try:

                SeqIO.write(

                    [
                        SeqRecord(

                            seq=Seq(record.seq.replace('\n', '').replace("\r", "")).replace("\t", ""),
                            id=record.id,
                            name=str(record.name).replace(" ", ""),
                            description=record.description,
                            dbxrefs=[record.dbxrefs],
                            annotations={
                                'molecule_type': obj.molecule,
                            },
                            letter_annotations={
                                'phred_quality': [record.letter_annotations["phred_quality"]]
                            }
                        )

                        for record in records

                    ],

                    out_handle.file.path,

                    obj.writing_file_format
                )

            except:

                try:

                    SeqIO.convert(
                        in_file=in_handle.file.path,
                        in_format=obj.reading_file_format,
                        out_file=out_handle.file.path,
                        out_format=obj.writing_file_format,
                        molecule_type=obj.molecule
                    )

                except ValueError as err:
                    messages.error(request, str(err))
                    return HttpResponseRedirect(request.get_full_path())

            return HttpResponseRedirect(reverse("bioinformatic:download"))

    return render(request, "bioinformatic/form.html", {'form': form, 'title': _("DOSYA DÖNÜŞTÜRME")})
