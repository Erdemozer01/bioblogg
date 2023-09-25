import os

from django.db.models import Q
from django.shortcuts import *
from django.views import generic
from django.contrib import messages
from bioinformatic.forms import FileReadingForm, TranslateForm, AlignmentForm
from bioinformatic.models import BioinformaticModel, RecordModel
from Bio import SeqIO
from Bio.SeqUtils import GC, gc_fraction
import pandas as pd
import plotly.express as px
from django_plotly_dash import DjangoDash
from dash import html, dcc, Input, Output, callback
from pathlib import Path
from Bio.Seq import Seq
from Bio import Align
from Bio.Align import substitution_matrices

BASE_DIR = Path(__file__).resolve().parent.parent.parent


def alignment_score(request):
    external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

    alignment_score = DjangoDash('alignment_score', external_stylesheets=external_stylesheets)

    alignment_score.layout = html.Div(

        [

            html.H4('ALİGNMENT SKOR'),

            html.A('BİYOİNFORMATİK ANASAYFA', href=HttpResponseRedirect(reverse("bioinformatic:home")).url,
                   style={'float': 'right'}),

            html.Div(
                [
                    html.P("Alignment Mod"),

                    dcc.Dropdown(id="mod", options={
                        'local': 'LOCAL',
                        'global': 'GLOBAL'
                    }, searchable=True),

                    html.P("Alignment MATRİX", style={'marginTop': '10px'}),

                    dcc.Dropdown(id="mat", options={
                        'BENNER22': 'BENNER22',
                        'BENNER6': 'BENNER6',
                        'BENNER74': 'BENNER74',
                        'BLOSUM45': 'BLOSUM45',
                        'BLOSUM50': 'BLOSUM50',
                        'BLOSUM62': 'BLOSUM62',
                        'BLOSUM80': 'BLOSUM80',
                        'BLOSUM90': 'BLOSUM90',
                        'DAYHOFF': 'DAYHOFF',
                        'FENG': 'FENG',
                        'HOXD70': 'HOXD70',
                        'JOHNSON': 'JOHNSON',
                        'JONES': 'JONES',
                        'LEVIN': 'LEVIN',
                        'MCLACHLAN': 'MCLACHLAN',
                        'MDM78': 'MDM78',
                        'PAM250': 'PAM250',
                        'PAM30': 'PAM30',
                        'PAM70': 'PAM70',
                        'RAO': 'RAO',
                        'RISLER': 'RISLER',
                        'SCHNEIDER': 'SCHNEIDER',
                        'STR': 'STR',
                        'TRANS': 'TRANS',
                    }, searchable=True),

                    html.P("Hedef", style={'marginTop': '10px'}),

                    dcc.Textarea(id="seq1", placeholder="Sekans Giriniz",
                                 style={'marginRight': '10px', 'width': '100%', 'height': '100px'}),

                    html.P("Sorgu", style={'marginTop': '10px'}),

                    dcc.Textarea(id="seq2", placeholder="Sekans Giriniz",
                                 style={'marginRight': '10px', 'width': '100%', 'height': '100px'}),

                ]
            ),

            html.Div(id="output"),

        ], style={'marginTop': "2%"}
    )

    @alignment_score.callback(
        Output("output", "children"),
        Input("mod", "value"),
        Input("mat", "value"),
        Input("seq1", "value"),
        Input("seq2", "value"),
    )
    def update_output(mod, mat, seq1, seq2):
        seq1 = str(seq1).replace(" ", '')
        seq1 = str(seq1).replace("\n", "")
        seq1 = str(seq1).replace("\t", "")
        seq1 = str(seq1).replace("\r", "")
        seq1 = str(seq1).replace("0", "")
        seq1 = str(seq1).replace("1", "")
        seq1 = str(seq1).replace("2", "")
        seq1 = str(seq1).replace("3", "")
        seq1 = str(seq1).replace("4", "")
        seq1 = str(seq1).replace("5", "")
        seq1 = str(seq1).replace("6", "")
        seq1 = str(seq1).replace("7", "")
        seq1 = str(seq1).replace("8", "")
        seq1 = str(seq1).replace("9", "")
        seq1 = str(seq1).upper()

        seq2 = str(seq2).replace(" ", "")
        seq2 = str(seq2).replace("\n", "")
        seq2 = str(seq2).replace("\t", "")
        seq2 = str(seq2).replace("\r", "")
        seq2 = str(seq2).replace("0", "")
        seq2 = str(seq2).replace("1", "")
        seq2 = str(seq2).replace("2", "")
        seq2 = str(seq2).replace("3", "")
        seq2 = str(seq2).replace("4", "")
        seq2 = str(seq2).replace("5", "")
        seq2 = str(seq2).replace("6", "")
        seq2 = str(seq2).replace("7", "")
        seq2 = str(seq2).replace("8", "")
        seq2 = str(seq2).replace("9", "")
        seq2 = str(seq2).upper()

        aligner = Align.PairwiseAligner()

        aligner.mode = str(mod)

        aligner.substitution_matrix = substitution_matrices.load(str(mat))

        alignments = aligner.align(seq1, seq2)

        alignment = alignments[0]

        count = alignment.counts()

        values = f'Boşluk : {count.gaps}  Benzerlik : {count.identities},  Eşleşmeyen: {count.mismatches},\n\nMatriks: {mat}, \n\n {alignment.substitutions}\n\nMOD: {str(mod).upper()} \n\nScore: {alignments.score}\n\nAlignment:\n\n{alignment}'

        return html.Div(
            [

                html.P(
                    f"Hedef sekans uzunluğu: {len(seq1)}, %GC: {GC(seq1)}, Sorgu sekans uzunluğu: {len(seq2)}, %GC: {GC(seq2)}",
                    style={'marginTop': '20px'}
                ),

                html.Hr(),

                html.Label("Alignments",  style={'marginBottom': '-35px'}),

                html.Button("İNDİR", id="btn-download-align",
                            style={'float': 'right', 'marginBottom': '20px'},
                            className='button-primary'),
                dcc.Download(id="download-align"),

                html.Div(
                    [
                        dcc.Textarea(id='align-textarea', value=values,
                                     style={'width': '100%', 'height': '500px', 'position': 'relative'}),
                    ]
                ),



                html.Hr(),
                html.Br(),

            ], id='output',

        )

    @alignment_score.callback(
        Output("download-align", "data"),
        Input("btn-download-align", "n_clicks"),
        Input('align-textarea', 'value'),
        prevent_initial_call=True,
    )
    def func(n_clicks, values):
        return dict(content=values, filename="alignment.txt")

    return HttpResponseRedirect("/laboratory/bioinformatic/app/alignment_score/")


def stats_view(request, user):
    global context
    if request.user.is_anonymous:
        from django.conf import settings
        messages.error(request, "Lütfen Giriş Yapınız")
        return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))

    object_list = RecordModel.objects.filter(records__user=request.user, records__tool="DOSYA OKUMA")

    if object_list.exists():

        df = pd.DataFrame(
            {
                "name": [dna.name for dna in object_list],
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
            to_stop = [i.to_stop]
        protein = [i.pro_seq_len for i in self.get_queryset()]
        context["title"] = "Sonuçlar"
        context["count"] = RecordModel.objects.filter(records__user=self.request.user,
                                                      records__tool="DOSYA OKUMA").count()
        context["to_stop"] = to_stop[0]

        if None in protein:
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

                for record in file_read:

                    if record.features:

                        for feature in record.features:

                            if feature.type == "CDS":

                                if feature.qualifiers.get('translation') is None:

                                    obj.record_content.create(
                                        record_id=record.id,
                                        name=record.name,
                                        sequence=record.seq,
                                        seq_len=len(record.seq),
                                        gene=feature.qualifiers.get('gene'),
                                        db_xrefs=feature.qualifiers.get('db_xref'),
                                        taxonomy=record.annotations['taxonomy'],
                                        description=record.description,
                                        gc=GC(record.seq).__round__(2),
                                    )

                                else:

                                    obj.record_content.create(
                                        record_id=record.id,
                                        name=record.name,
                                        taxonomy=record.annotations['taxonomy'],
                                        description=record.description,
                                        gene=feature.qualifiers.get('gene'),
                                        db_xrefs=feature.qualifiers.get('db_xref'),
                                        protein_sequence=feature.qualifiers.get('translation')[0],
                                        pro_seq_len=len(feature.qualifiers.get('translation')[0]),
                                        sequence=record.seq,
                                        seq_len=len(record.seq),
                                        gc=GC(record.seq).__round__(2),
                                    )

                            else:

                                obj.record_content.create(
                                    record_id=record.id,
                                    name=record.name,
                                    sequence=record.seq,
                                    seq_len=len(record.seq),
                                    db_xrefs=record.dbxrefs,
                                    taxonomy=record.annotations['taxonomy'],
                                    description=record.description,
                                    gc=GC(record.seq).__round__(2),
                                )

                    else:

                        for record in file_read:
                            obj.record_content.create(
                                record_id=record.id,
                                name=record.name,
                                description=record.description,
                                sequence=record.seq,
                                db_xrefs=record.dbxrefs,
                                annotations=record.annotations,
                                features=record.features,
                                gc=GC(record.seq).__round__(2),
                                seq_len=len(record.seq),
                            )

            except UnicodeDecodeError:
                BioinformaticModel.objects.filter(user=request.user).delete()
                return render(request, "exception/page-404.html", {'msg': 'Hatalı dosya türü', 'url': request.path})
            file_obj.delete()
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
    if request.user.is_anonymous:
        from django.conf import settings
        messages.error(request, "Lütfen Giriş Yapınız")
        return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))

    form = TranslateForm(request.POST or None)
    object_list = BioinformaticModel.objects.filter(user=request.user, tool="DOSYA OKUMA")
    record_obj = RecordModel.objects.filter(records__user=request.user, records__tool="DOSYA OKUMA")

    if request.method == "POST":

        if form.is_valid():

            trans_table = form.cleaned_data["trans_table"]
            to_stop = form.cleaned_data["to_stop"]

            if to_stop is True:

                for object in object_list:
                    object.to_stop = to_stop
                    object.trans_table = trans_table
                    object.save()

                for record in record_obj:
                    record.protein_sequence = Seq(record.sequence).translate(table=trans_table)
                    record.pro_seq_len = len(Seq(record.sequence).translate(table=trans_table))
                    record.save()

            else:
                for object in object_list:
                    object.to_stop = to_stop
                    object.trans_table = trans_table
                    object.save()

                for record in record_obj:
                    record.protein_sequence = Seq(record.sequence).translate(table=trans_table).replace("*", "")
                    record.pro_seq_len = len(Seq(record.sequence).translate(table=trans_table).replace("*", ""))
                    record.save()

            return HttpResponseRedirect(
                reverse("bioinformatic:file_reading_results", kwargs={'user': request.user}))

    return render(request, "bioinformatic/form.html", {"form": form, "title": "Proteinleri Topla"})
