import os

from django.db.models import Q
from django.shortcuts import *
from django.views import generic
from django.contrib import messages
from bioinformatic.forms import FileReadingForm, TranslateForm, AlignmentForm
from bioinformatic.models import BioinformaticModel, RecordModel
from Bio import SeqIO
from Bio.SeqUtils import GC
import pandas as pd
import plotly.express as px
from django_plotly_dash import DjangoDash
from dash import html, dcc, Input, Output, callback
from pathlib import Path
from Bio.Seq import Seq
from Bio import Align
from Bio.Align import substitution_matrices

BASE_DIR = Path(__file__).resolve().parent.parent.parent


def alignment_score(request, user):

    if request.user.is_anonymous:
        from django.conf import settings
        messages.error(request, "Lütfen Giriş Yapınız")
        return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))

    form = AlignmentForm(request.POST or None)

    if request.method == "POST":

        if form.is_valid():

            try:
                aligner = Align.PairwiseAligner()

                seq1 = form.cleaned_data['seq1'].replace("\n", "").replace("\r", "")

                seq2 = form.cleaned_data['seq2'].replace("\n", "").replace("\r", "")

                aligner.mode = form.cleaned_data['mode']

                aligner.substitution_matrix = substitution_matrices.load(form.cleaned_data['matrix'])

                mod = form.cleaned_data['mode']

                matrix = form.cleaned_data['matrix']

                alignments = aligner.align(seq1, seq2)

                alignment = alignments[0]

                count = alignment.counts()

                app = DjangoDash('alignment')

                app.layout = html.Div([
                    dcc.Textarea(
                        id='textarea-example',
                        value=f'\nBoşluk : {count.gaps}, Benzerlik : {count.identities}, Eşleşmeyen: {count.mismatches},\n\n {alignment.substitutions} \n\nMatriks: {matrix} \n\nMOD: {mod} \n\nScore: {alignment.score} \n\n' + "Alignment: \n\n" + f"{alignment}",
                        style={'width': '100%', 'height': 300},
                        className="form-control mb-2",
                    ),
                ])

            except ValueError:
                return render(request, "exception/page-404.html", {'msg': 'Sekansınızı kontrol ediniz.'})

            return render(request, "bioinformatic/reading/alignment.html")

        else:

            form = AlignmentForm()

    return render(request, "bioinformatic/form.html", {'form': form, 'title': 'Alignment'})


def stats_view(request, user):
    global context
    if request.user.is_anonymous:
        from django.conf import settings
        messages.error(request, "Lütfen Giriş Yapınız")
        return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))

    object_list = RecordModel.objects.filter(user=request.user, tool="DOSYA OKUMA")

    if object_list.exists():

        df = pd.DataFrame(
            {
                "name": [dna.molecule_id for dna in object_list],
                "DNA": [dna.seq_len for dna in object_list],
                "PROTEİN": [pro.pro_seq_len for pro in object_list],
                "%GC": [gc.gc for gc in object_list],
            }
        )

        scatter_map = DjangoDash(f"l{request.user}r")

        scatter_map.layout = html.Div([
            dcc.Graph(
                figure=px.scatter(
                    data_frame=df, x="DNA", y="PROTEİN",
                    title="DNA - PROTEİN Lineer Regrasyon", trendline="ols",
                    labels={
                        "dna_seq_len": "DNA SEKANS UZUNLUĞU",
                        "pro_seq_len": "PROTEİN SEKANS UZUNLUĞU",
                        'name': 'Canlı'
                    },
                )),
        ])

        context = {

            'stat': df.describe().round(3).to_html(
                classes="table table-hover shadow-soft rounded-lg text-center",
                border=False, header="İstatistik", justify="center"
            ),

            'table_cov': df.cov().round(3).to_html(
                classes="table table-hover shadow-soft rounded-lg text-center",
                border=False, header="Kovaryans", justify="center"
            ),

            'table_cor': df.corr().round(3).to_html(
                classes="table table-hover shadow-soft rounded-lg text-center",
                border=False, header="Korelasyon", justify="center"

            ),

        }

    else:

        return render(request, "exception/page-404.html", {'msg': "Verilere Ulaşılamadı"})

    return render(request, "bioinformatic/stats/stats.html", context)


class FileReadingResultView(generic.ListView):
    template_name = "bioinformatic/reading/result.html"
    model = RecordModel
    paginate_by = 10

    def get(self, request, *args, **kwargs):
        if request.user.is_anonymous:
            from django.conf import settings
            messages.error(request, "Lütfen Giriş Yapınız")
            return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))
        return super().get(request, *args, **kwargs)

    def get_queryset(self):
        search = self.request.GET.get('search', False)
        if search:
            return RecordModel.objects.filter(
                Q(description__icontains=search, records__user=self.request.user, records__tool="DOSYA OKUMA"))
        else:
            return RecordModel.objects.filter(records__user=self.request.user, records__tool="DOSYA OKUMA")

    def get_context_data(self, *, object_list=None, **kwargs):
        global file_format, protein, to_stop
        context = super().get_context_data(**kwargs)
        for i in BioinformaticModel.objects.filter(user=self.request.user, tool="DOSYA OKUMA"):
            file_format = [i.reading_file_format]
            protein = [j.protein_sequence for j in
                       i.record_content.filter(records__user=self.request.user, records__tool="DOSYA OKUMA")]
            to_stop = [i.to_stop]

        context["title"] = "Sonuçlar"
        context["count"] = RecordModel.objects.filter(records__user=self.request.user, records__tool="DOSYA OKUMA").count()

        if protein[0] == '':
            context["protein"] = "protein yok"
        else:
            context["protein"] = "protein var"

        page = context['page_obj']
        paginator = page.paginator
        pagelist = paginator.get_elided_page_range(page.number, on_each_side=2, on_ends=2)

        context['pagelist'] = pagelist

        return context


def file_reading(request, user):
    global file_read

    if request.user.is_anonymous:
        from django.conf import settings
        messages.error(request, "Lütfen Giriş Yapınız")
        return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))

    if BioinformaticModel.objects.filter(user=request.user, tool="DOSYA OKUMA").exists():
        BioinformaticModel.objects.filter(user=request.user, tool="DOSYA OKUMA").delete()

    form = FileReadingForm(request.POST or None, request.FILES or None)

    if request.method == "POST":
        if form.is_valid():

            file_format = form.cleaned_data["reading_file_format"]
            read_file = form.cleaned_data["file"]
            molecule = form.cleaned_data["molecule"]

            obj = BioinformaticModel.objects.create(
                user=request.user,
                molecule=molecule,
                reading_file_format=file_format,
                tool="DOSYA OKUMA"
            )

            file_obj = obj.records_files.create(
                file=read_file
            )

            if read_file.size > 5242880:
                obj.delete()
                messages.error(request, "Dosya boyutu 5mb dan büyüktür.")
                return redirect(request.META['HTTP_REFERER'])

            file_read = SeqIO.parse(file_obj.file.path, file_format)

            try:

                if file_format == "gb" or "genbank" or "genbank-cds":

                    for record in file_read:

                        for feature in record.features:

                            if feature.type == "CDS":

                                if feature.qualifiers.get('translation') is not None:

                                    obj.records.create(
                                        record_id=record.id,
                                        name=record.name,
                                        sequence=record.seq,
                                        seq_len=len(record.seq),
                                        organism=feature.qualifiers.get('organism'),
                                        gene=feature.qualifiers.get('gene'),
                                        dbxrefs=feature.qualifiers.get('db_xref'),
                                        source=feature.qualifiers.get('source'),
                                        taxonomy=record.annotations['taxonomy'],
                                        keywords=record.annotations['keywords'],
                                        description=record.description,
                                        gc=GC(record.seq).__round__(2),
                                    )

                                else:

                                    obj.records.create(
                                        organism=feature.qualifiers.get('organism'),
                                        gene=feature.qualifiers.get('gene'),
                                        dbxrefs=feature.qualifiers.get('db_xref'),
                                        source=feature.qualifiers.get('source'),
                                        protein_sequence=feature.qualifiers.get('translation')[0],
                                        protein_id=feature.qualifiers.get('protein_id')[0],
                                        pro_seq_len=len(feature.qualifiers.get('translation')[0]),
                                        description=record.description,
                                        sequence=record.seq,
                                        seq_len=len(record.seq),
                                        gc=GC(record.seq).__round__(2),
                                    )

                            else:

                                obj.records.create(
                                    record_id=record.id,
                                    name=record.name,
                                    sequence=record.seq,
                                    seq_len=len(record.seq),
                                    taxonomy=record.annotations['taxonomy'],
                                    keywords=record.annotations['keywords'],
                                    description=record.description,
                                    gc=GC(record.seq).__round__(2),
                                )

                else:

                    for record in file_read:
                        obj.records.create(
                            record_id=record.id,
                            name=record.name,
                            description=record.description,
                            sequence=record.seq,
                            dbxrefs=record.dbxrefs,
                            annotations=record.annotations,
                            features=record.features,
                            gc=GC(record.seq).__round__(2),
                            seq_len=len(record.seq),
                        )

            except UnicodeDecodeError:
                BioinformaticModel.objects.get(user=request.user, tool="").delete()
                return render(request, "exception/page-404.html", {'msg': 'Hatalı dosya türü', 'url': request.path})

            return HttpResponseRedirect(reverse("bioinformatic:file_reading_results", kwargs={'user': request.user}))

        else:

            form = FileReadingForm()

    return render(request, "bioinformatic/form.html", {'form': form, 'title': 'DOSYA OKUMASI'})


class FileReadDetailView(generic.DetailView):
    template_name = "bioinformatic/reading/detail.html"
    model = RecordModel

    def get_queryset(self):
        return RecordModel.objects.filter(records__user=self.request.user, records__tool="DOSYA OKUMA")

    def get(self, request, *args, **kwargs):
        if request.user.is_anonymous:
            from django.conf import settings
            messages.error(request, "Lütfen Giriş Yapınız")
            return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))
        return super().get(request, *args, **kwargs)

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['A'] = self.object.sequence.count("A")
        context['T'] = self.object.sequence.count("T")
        context['G'] = self.object.sequence.count("G")
        context['C'] = self.object.sequence.count("C")
        return context


def ProteinPickView(request, user):
    form = TranslateForm(request.POST or None)
    object_list = BioinformaticModel.objects.filter(user=request.user, tool="okuma")

    if request.method == "POST":

        if form.is_valid():

            trans_table = form.cleaned_data["trans_table"]
            to_stop = form.cleaned_data["to_stop"]

            if to_stop is True:
                for object in object_list:
                    object.to_stop = to_stop
                    object.pro_seq = Seq(object.seq).translate(table=trans_table)
                    object.pro_seq_len = len(Seq(object.seq).translate(table=trans_table))
                    object.save()

            else:
                for object in object_list:
                    object.to_stop = to_stop
                    object.pro_seq = Seq(object.seq).translate(table=trans_table).replace("*", "")
                    object.pro_seq_len = len(Seq(object.seq).translate(table=trans_table).replace("*", ""))
                    object.save()

            return HttpResponseRedirect(
                reverse("bioinformatic:file_reading_results", kwargs={'user': request.user}))

    return render(request, "bioinformatic/form.html", {"form": form, "title": "Proteinleri Topla"})
