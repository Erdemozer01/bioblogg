import os

from django.db.models import Q
from django.shortcuts import *
from django.views import generic
from django.contrib import messages
from bioinformatic.forms import FileReadingForm, TranslateForm, FileResulSelect, AlignmentForm
from bioinformatic.models import BioinformaticAnalizModel
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

    object_list = BioinformaticAnalizModel.objects.filter(user=request.user, tool="okuma")

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
    model = BioinformaticAnalizModel
    paginate_by = 10
    form_class = FileResulSelect

    def get(self, request, *args, **kwargs):
        if request.user.is_anonymous:
            from django.conf import settings
            messages.error(request, "Lütfen Giriş Yapınız")
            return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))
        return super().get(request, *args, **kwargs)

    def get_queryset(self):
        search = self.request.GET.get('search', False)
        if search:
            return BioinformaticAnalizModel.objects.filter(
                Q(description__icontains=search, user=self.request.user, tool="okuma"))
        else:
            return BioinformaticAnalizModel.objects.filter(user=self.request.user, tool="okuma")

    def get_context_data(self, *, object_list=None, **kwargs):
        global file_format, protein, to_stop
        context = super().get_context_data(**kwargs)
        for i in BioinformaticAnalizModel.objects.filter(user=self.request.user, tool="okuma"):
            file_format = [i.reading_file_format]
            protein = [i.pro_seq]
            to_stop = [i.to_stop]

        context["title"] = "Sonuçlar"
        context["count"] = BioinformaticAnalizModel.objects.filter(user=self.request.user, tool="okuma").count()

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

    if BioinformaticAnalizModel.objects.filter(user=request.user).exists():
        BioinformaticAnalizModel.objects.filter(user=request.user).delete()

    form = FileReadingForm(request.POST or None, request.FILES or None)

    if request.method == "POST":
        if form.is_valid():

            file_format = form.cleaned_data["reading_file_format"]
            read_file = form.cleaned_data["read_file"]
            molecule = form.cleaned_data["molecule"]

            obj = BioinformaticAnalizModel.objects.create(
                user=request.user,
                read_file=read_file,
                molecule=molecule,
                file_format=file_format
            )

            if read_file.size > 5242880:
                obj.delete()
                messages.error(request, "Dosya boyutu 5mb dan büyüktür.")
                return redirect(request.META['HTTP_REFERER'])

            file_read = SeqIO.parse(obj.read_file.path, file_format)

            try:

                if molecule == "DNA":

                    for file_content in file_read:
                        BioinformaticAnalizModel.objects.create(
                            user=request.user,
                            file_format=file_format,
                            tool="okuma",
                            molecule=molecule,
                            molecule_id=file_content.id,
                            seq=file_content.seq,
                            annotations=file_content.annotations,
                            description=file_content.description,
                            features=file_content.features,
                            gc=GC(file_content.seq).__round__(2),
                            seq_len=len(file_content.seq),
                            organism=file_content.name
                        )

                elif molecule == "protein":

                    for record in file_read:

                        for feature in record.features:

                            print(feature)

                            if feature.type == "CDS":

                                if feature.qualifiers.get('translation') is not None:

                                    print(feature.qualifiers.get('organism'))

                                    BioinformaticAnalizModel.objects.create(
                                        user=request.user,
                                        tool="okuma",
                                        file_format=file_format,
                                        molecule=molecule,
                                        seq=record.seq,
                                        seq_len=len(record.seq),
                                        pro_seq=feature.qualifiers.get('translation')[0],
                                        molecule_id=feature.qualifiers.get('protein_id')[0],
                                        pro_seq_len=len(feature.qualifiers.get('translation')[0]),
                                        organism=feature.qualifiers.get('organism'),
                                        taxonomy=record.annotations['taxonomy'],
                                        description=record.description,
                                        gc=GC(record.seq).__round__(2),
                                    )

                                else:

                                    BioinformaticAnalizModel.objects.create(
                                        user=request.user,
                                        tool="okuma",
                                        molecule=molecule,
                                        file_format=file_format,
                                        organism=record.annotations['organism'],
                                        taxonomy=record.annotations['taxonomy'],
                                        description=record.description,
                                        seq=record.seq,
                                        seq_len=len(record.seq),
                                        gc=GC(record.seq).__round__(2),
                                    )

                            else:

                                BioinformaticAnalizModel.objects.create(
                                    user=request.user,
                                    tool="okuma",
                                    molecule=molecule,
                                    file_format=file_format,
                                    organism=record.annotations['organism'],
                                    taxonomy=record.annotations['taxonomy'],
                                    description=record.description,
                                    seq=record.seq,
                                    seq_len=len(record.seq),
                                    gc=GC(record.seq).__round__(2),
                                )

            except UnicodeDecodeError:
                BioinformaticAnalizModel.objects.get(user=request.user, tool="").delete()
                return render(request, "exception/page-404.html", {'msg': 'Hatalı dosya türü', 'url': request.path})

            BioinformaticAnalizModel.objects.get(user=request.user, tool="").delete()

            return HttpResponseRedirect(reverse("bioinformatic:file_reading_results", kwargs={'user': request.user}))

        else:

            form = FileReadingForm()

    return render(request, "bioinformatic/form.html", {'form': form, 'title': 'DOSYA OKUMASI'})


class FileReadDetailView(generic.DetailView):
    template_name = "bioinformatic/reading/detail.html"
    model = BioinformaticAnalizModel

    def get_queryset(self):
        return BioinformaticAnalizModel.objects.filter(user=self.request.user, tool="okuma")

    def get(self, request, *args, **kwargs):
        if request.user.is_anonymous:
            from django.conf import settings
            messages.error(request, "Lütfen Giriş Yapınız")
            return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))
        return super().get(request, *args, **kwargs)

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['A'] = self.object.seq.count("A")
        context['T'] = self.object.seq.count("T")
        context['G'] = self.object.seq.count("G")
        context['C'] = self.object.seq.count("C")
        return context


def ProteinPickView(request, user):
    form = TranslateForm(request.POST or None)
    object_list = BioinformaticAnalizModel.objects.filter(user=request.user, tool="okuma")
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
