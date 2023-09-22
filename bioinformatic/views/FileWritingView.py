import os.path
from django.views import generic
from django.shortcuts import *
from pathlib import Path
from django.contrib import messages
from bioinformatic.forms.writing import SelectWritingForm, FileWritingForm
from bioinformatic.models.bioinformatic import BioinformaticModel, RecordModel
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

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


class FileWritingListView(generic.ListView):
    template_name = "bioinformatic/writing/list.html"
    model = RecordModel
    paginate_by = 10

    def get(self, request, *args, **kwargs):
        if request.user.is_anonymous:
            from django.conf import settings
            messages.error(request, "Lütfen Giriş Yapınız")
            return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))
        return super().get(request, *args, **kwargs)

    def get_queryset(self):
        return RecordModel.objects.filter(records__user=self.request.user, records__tool="DOSYA YAZMA")

    def get_context_data(self, *, object_list=None, **kwargs):
        context = super().get_context_data(**kwargs)
        context['count'] = self.object_list.count()
        context['format'] = BioinformaticModel.objects.filter(user=self.request.user,
                                                              tool="DOSYA YAZMA").last().writing_file_format
        context['title'] = f"{context['format'].upper()} KAYITLARI"
        return context


def FileWritingView(request, format, user):
    if request.user.is_anonymous:
        from django.conf import settings
        messages.error(request, "Lütfen Giriş Yapınız")
        return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))

    form = FileWritingForm(request.POST or None)

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

            object_list = BioinformaticModel.objects.filter(user=request.user, tool="DOSYA YAZMA")

            for obj in object_list:
                obj.record_content.create(
                    molecule=form.cleaned_data['molecule'],
                    molecule_id=form.cleaned_data['molecule_id'].replace("|", ""),
                    name=form.cleaned_data['name'],
                    description=form.cleaned_data['description'],
                    db_xrefs=form.cleaned_data['db_xrefs'],
                    annotations=form.cleaned_data['annotations'],
                    source=form.cleaned_data['source'],
                    organism=form.cleaned_data['organism'],
                    keywords=form.cleaned_data['keywords'],
                    accession=form.cleaned_data['accessions'],
                    sequence=str(sekans)
                )

            return HttpResponseRedirect(
                reverse("bioinformatic:file_writing_list", kwargs={'format': format, 'user': request.user}))

        else:

            form = FileWritingForm()

    return render(
        request, 'bioinformatic/writing/form.html',
        {
            'form': form,
            'title': f'{format.upper()} DOSYASI OLUŞTURMA',
            'format': format
        }
    )


def CreateFileView(request, user, format):
    obj_list = RecordModel.objects.filter(
        records__user=request.user,
        records__tool='DOSYA YAZMA',
        records__writing_file_format=format
    )

    file_path = os.path.join(BASE_DIR, 'media', f'{user}_{format}_file.{format}')

    SeqIO.write(

        [
            SeqRecord(

                seq=Seq(record.sequence.replace('\n', '').replace("\r", "")).replace("\t", ""),
                id=record.molecule_id,
                name=str(record.name).replace(" ", ""),
                description=record.description,
                dbxrefs=[record.db_xrefs],
                annotations={
                    'molecule_type': [record.molecule],
                    'source': [record.source],
                    'keywords': [str(record.keywords)],
                    'organism': [record.organism],
                    'accession': [record.accession],
                    'MDAT': [record.records.created],
                },

            )

            for record in obj_list

        ], file_path, format
    )

    return HttpResponseRedirect(reverse('bioinformatic:download'))


class RecordDetailView(generic.DetailView):
    template_name = "bioinformatic/writing/detail.html"
    model = RecordModel


def RecordDeleteView(request, pk):
    RecordModel.objects.get(pk=pk).delete()
    return redirect(request.META['HTTP_REFERER'])
