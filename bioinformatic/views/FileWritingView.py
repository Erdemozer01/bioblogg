import os.path
import datetime
from pathlib import Path
from django.contrib import messages
from bioinformatic.forms.writing import SelectWritingForm, GenbankWritingForm, FastaWritingForm
from bioinformatic.models.bioinformatic import BioinformaticModel, RecordModel
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from django.shortcuts import *

BASE_DIR = Path(__file__).resolve().parent.parent.parent


def file_writing_format_select(request, user):
    global file_format

    if request.user.is_anonymous:
        from django.conf import settings
        messages.error(request, "Lütfen Giriş Yapınız")
        return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))

    if BioinformaticModel.objects.filter(user=request.user, tool="DOSYA YAZMA").exists():
        BioinformaticModel.objects.filter(user=request.user, tool="DOSYA YAZMA").delete()

    form = SelectWritingForm(request.POST or None)

    if request.method == "POST":
        if form.is_valid():
            file_format = form.cleaned_data['writing_file_format']

            BioinformaticModel.objects.create(
                user=request.user,
                tool='DOSYA YAZMA',
                writing_file_format=file_format,
            )

            return HttpResponseRedirect(
                reverse("bioinformatic:file_writing", kwargs={'user': request.user, 'format': file_format}))

    return render(request, 'bioinformatic/form.html', {'form': form, 'title': 'DOSYA OLUŞTURMA'})


def FileWritingView(request, user, format):
    if request.user.is_anonymous:
        from django.conf import settings
        messages.error(request, "Lütfen Giriş Yapınız")
        return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))

    if format == 'fasta':
        object_list = BioinformaticModel.objects.filter(user=request.user, tool="DOSYA YAZMA")
        rec_obj = RecordModel.objects.filter(records__user=request.user, records__tool="DOSYA YAZMA",
                                             records__writing_file_format=format)
        form = FastaWritingForm(request.POST or None)

        if request.method == "POST":
            if form.is_valid():

                sekans = form.cleaned_data['sequence'].upper()

                sekans = sekans.replace(" ", "")
                sekans = sekans.replace("\n", "")
                sekans = sekans.replace("\t", "")
                sekans = sekans.replace("\r", "")
                sekans = sekans.replace("0", "")
                sekans = sekans.replace("1", "")
                sekans = sekans.replace("2", "")
                sekans = sekans.replace("3", "")
                sekans = sekans.replace("4", "")
                sekans = sekans.replace("5", "")
                sekans = sekans.replace("6", "")
                sekans = sekans.replace("7", "")
                sekans = sekans.replace("8", "")
                sekans = sekans.replace("9", "")

                for obj in object_list:
                    obj.record_content.create(
                        molecule=form.cleaned_data['molecule'],
                        molecule_id=form.cleaned_data['molecule_id'].replace("|", ""),
                        name=form.cleaned_data['name'],
                        description=form.cleaned_data['description'],
                        sequence=str(sekans)
                    )

                return HttpResponseRedirect(request.META.get('HTTP_REFERER'))

        return render(request, "bioinformatic/writing/list.html",
                      {'form': form, 'rec_obj': rec_obj, 'format': format,
                       'title': f'{format} dosyası oluşturma'.upper()})

    elif format == 'genbank':
        object_list = BioinformaticModel.objects.filter(user=request.user, tool="DOSYA YAZMA")
        rec_obj = RecordModel.objects.filter(records__user=request.user, records__tool="DOSYA YAZMA",
                                             records__writing_file_format=format)
        form = GenbankWritingForm(request.POST or None)

        if request.method == "POST":
            if form.is_valid():

                sekans = form.cleaned_data['sequence'].upper()

                sekans = sekans.replace(" ", "")
                sekans = sekans.replace("\n", "")
                sekans = sekans.replace("\t", "")
                sekans = sekans.replace("\r", "")
                sekans = sekans.replace("0", "")
                sekans = sekans.replace("1", "")
                sekans = sekans.replace("2", "")
                sekans = sekans.replace("3", "")
                sekans = sekans.replace("4", "")
                sekans = sekans.replace("5", "")
                sekans = sekans.replace("6", "")
                sekans = sekans.replace("7", "")
                sekans = sekans.replace("8", "")
                sekans = sekans.replace("9", "")

                for obj in object_list:
                    obj.record_content.create(
                        molecule=form.cleaned_data['molecule'],
                        molecule_id=form.cleaned_data['molecule_id'].replace("|", ""),
                        name=form.cleaned_data['name'],
                        description=form.cleaned_data['description'],
                        data_file_division='CON',
                        date=str(datetime.datetime.now().strftime("%d-%b-%Y")).upper(),
                        source=form.cleaned_data['source'],
                        taxonomy=form.cleaned_data['taxonomy'],
                        organism=form.cleaned_data['organism'],
                        keywords=".",
                        accession=form.cleaned_data['accessions'],
                        topology=form.cleaned_data["topology"],
                        sequence=str(sekans)
                    )

                return HttpResponseRedirect(request.META.get('HTTP_REFERER'))

        return render(request, "bioinformatic/writing/list.html",
                      {'form': form, 'rec_obj': rec_obj, 'format': format,
                       'title': f'{format} dosyası oluşturma'.upper()})


def CreateFileView(request, user, format):
    if request.user.is_anonymous:
        from django.conf import settings
        messages.error(request, "Lütfen Giriş Yapınız")
        return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))

    obj_list = RecordModel.objects.filter(
        records__user=request.user,
        records__tool='DOSYA YAZMA',
        records__writing_file_format=format
    )

    write_file_path = os.path.join(BASE_DIR, 'media', 'laboratory', f'{user}_{format}_file.{format}')

    SeqIO.write(

        [
            SeqRecord(

                seq=Seq(str(record.sequence)),
                id=record.molecule_id,
                name=str(record.name).replace(" ", ""),
                description=record.description,
                dbxrefs=record.db_xrefs,
                annotations={
                    'molecule_type': record.molecule,
                    'taxonomy': [record.taxonomy],
                    'data_file_division': 'CON',
                    'date': str(datetime.datetime.now().strftime("%d-%b-%Y")).upper(),
                    'source': record.source,
                    'keywords': [str(record.keywords)],
                    'organism': record.organism,
                    'accession': record.accession,
                    'topology': record.topology
                },

            )

            for record in obj_list

        ], write_file_path, format
    )

    return HttpResponseRedirect(reverse('bioinformatic:download_file'))


def RecordDeleteView(request, pk):
    RecordModel.objects.get(pk=pk).delete()
    return redirect(request.META['HTTP_REFERER'])
