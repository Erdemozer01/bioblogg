import os.path
from django.views import generic
from django.shortcuts import *
from pathlib import Path
from django.contrib import messages
from bioinformatic.forms.writing import SelectWritingForm, FileWritingForm
from bioinformatic.models.bioinformatic import BioinformaticAnalizModel
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

    if BioinformaticAnalizModel.objects.filter(user=request.user, tool="DOSYA YAZMA").exists():
        BioinformaticAnalizModel.objects.filter(user=request.user, tool="DOSYA YAZMA").delete()

    form = SelectWritingForm(request.POST or None)

    if request.method == "POST":
        if form.is_valid():

            file_format = form.cleaned_data['writing_file_format']
            if file_format == "fasta":

                BioinformaticAnalizModel.objects.create(
                    user=request.user,
                    tool='DOSYA YAZMA',
                    writing_file_format=file_format,
                )

            else:
                BioinformaticAnalizModel.objects.create(
                    user=request.user,
                    tool='DOSYA YAZMA',
                    writing_file_format=file_format,
                )

            return HttpResponseRedirect(
                reverse("bioinformatic:file_writing", kwargs={'user': request.user, 'format': file_format}))

    return render(request, 'bioinformatic/form.html', {'form': form, 'title': 'DOSYA OLUŞTURMA'})


class FileWritingListView(generic.ListView):
    template_name = "bioinformatic/writing/list.html"
    model = BioinformaticAnalizModel
    paginate_by = 10

    def get(self, request, *args, **kwargs):
        if request.user.is_anonymous:
            from django.conf import settings
            messages.error(request, "Lütfen Giriş Yapınız")
            return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))
        return super().get(request, *args, **kwargs)

    def get_queryset(self):
        return BioinformaticAnalizModel.objects.filter(user=self.request.user, tool="DOSYA YAZMA")[1:]

    def get_context_data(self, *, object_list=None, **kwargs):
        context = super().get_context_data(**kwargs)
        context['count'] = self.object_list.count()
        context['format'] = BioinformaticAnalizModel.objects.filter(user=self.request.user,
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

            BioinformaticAnalizModel.objects.create(
                user=request.user,
                tool="DOSYA YAZMA",
                writing_file_format=format,
                name=form.cleaned_data['name'],
                molecule=form.cleaned_data['molecule'],
                molecule_id=form.cleaned_data['molecule_id'].replace("|", ""),
                description=form.cleaned_data['description'],
                dbxrefs=form.cleaned_data['dbxrefs'],
                source=form.cleaned_data['source'],
                keywords=form.cleaned_data['keywords'],
                organism=form.cleaned_data['organism'],
                accession=form.cleaned_data['accession'],
                seq=str(
                    form.cleaned_data['seq'].upper().replace("\n", "").replace('\r', '').replace(" ", "").replace("\t",
                                                                                 "")),
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
    obj_list = BioinformaticAnalizModel.objects.filter(
        user=request.user,
        tool='DOSYA YAZMA',
        writing_file_format=format
    )[1:]

    file_path = os.path.join(BASE_DIR, 'media', f'{request.user}_{format}_file.{format}')

    SeqIO.write(

        [
            SeqRecord(

                seq=Seq(record.seq.replace('\n', '').replace("\r", "")).replace("\t", ""),
                id=record.molecule_id,
                name=str(record.name).replace(" ", ""),
                description=record.description,
                dbxrefs=[record.dbxrefs],
                annotations={
                    'molecule_type': [record.molecule],
                    'source': [record.source],
                    'keywords': [str(record.keywords)],
                    'organism': [record.organism],
                    'accession': [record.accession],
                    'MDAT': [record.created],
                },

            )

            for record in obj_list

        ], file_path, format
    )

    return HttpResponseRedirect(reverse('bioinformatic:download'))


class RecordDetailView(generic.DetailView):
    template_name = "bioinformatic/writing/detail.html"
    model = BioinformaticAnalizModel


def RecordDeleteView(request, pk):
    BioinformaticAnalizModel.objects.get(pk=pk).delete()
    return redirect(request.META['HTTP_REFERER'])
