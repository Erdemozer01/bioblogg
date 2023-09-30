import os
import gzip
from Bio import bgzf
from django.db.models import Q
from django.shortcuts import *
from django.views import generic
from django.contrib import messages
from bioinformatic.forms import FileReadingForm, TranslateForm, AlignmentForm
from bioinformatic.models import BioinformaticModel, RecordModel, FileModel
from Bio import SeqIO
from Bio.SeqUtils import GC, gc_fraction
import pandas as pd
import pandas as gen
import plotly.express as px
from django_plotly_dash import DjangoDash
from dash import Dash, html, dcc, dash_table
from pathlib import Path
from Bio.Seq import Seq
import dash_ag_grid as dag
from collections import defaultdict

BASE_DIR = Path(__file__).resolve().parent.parent.parent


def file_reading(request, user):
    if request.user.is_anonymous:
        from django.conf import settings
        messages.error(request, "Lütfen Giriş Yapınız")
        return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))

    external_stylesheets = ['https://cdn.jsdelivr.net/npm/bootstrap@5.3.2/dist/css/bootstrap.min.css']
    external_scripts = ['https://cdn.jsdelivr.net/npm/bootstrap@5.3.2/dist/js/bootstrap.bundle.min.js']

    if BioinformaticModel.objects.filter(user=request.user, tool="DOSYA OKUMA").exists():
        BioinformaticModel.objects.filter(user=request.user, tool="DOSYA OKUMA").delete()

    form = FileReadingForm(request.POST or None, request.FILES or None)

    if request.method == "POST":

        if request.FILES['file'].size > 100 * 1024 * 1024:
            messages.error(request, 'Dosya boyutu en en fazla 5 mb dan fazladır.')
            return HttpResponseRedirect(request.path)

        if form.is_valid():

            file_format = form.cleaned_data["reading_file_format"]
            file = form.cleaned_data["file"]
            molecule = form.cleaned_data["molecule"]

            obj = BioinformaticModel.objects.create(
                user=request.user,
                molecule=molecule,
                reading_file_format=file_format,
                tool="DOSYA OKUMA",
            )

            file_obj = obj.records_files.create(file=file)

            handle = open(file_obj.file.path, 'r', encoding='utf-8')

            file_reading_dash_app = DjangoDash(
                name=f"{file_format}-dosya-sonuc",
                external_stylesheets=external_stylesheets,
                external_scripts=external_scripts,
                title=f"{file_format} dosyası okuması sonuçları".upper()
            )

            try:
                records = SeqIO.parse(handle, file_format)

            except UnicodeDecodeError:
                try:
                    records = SeqIO.parse(gzip.open(file_obj.file.path, 'rt', encoding='utf-8'), file_format)

                except:
                    obj.delete()
                    return render(request, "exception/page-404.html", {'msg': "Hatalı dosya türü"})
            try:
                if file_format == "fasta":

                    df_TABLE = pd.DataFrame(
                        {
                            'id': str(rec.id),
                            'description': str(rec.description),
                            'seq': str(rec.seq),
                            'seq_len': len(str(rec.seq)),
                            'gc': gc_fraction(rec.seq) * 100
                        } for rec in records)

                    file_reading_dash_app.layout = html.Div(
                        children=[

                            html.H4("Fasta dosyası okuması", className="fw-bolder text-primary m-2 "),

                            html.Hr(className="border border-danger"),

                            dag.AgGrid(
                                id="df-table",
                                style={'width': '100%'},
                                rowData=df_TABLE.to_dict("records"),
                                columnDefs=[
                                    {'field': 'id', 'headerName': 'İD', 'filter': True},
                                    {'field': 'description', 'headerName': 'Tanım', 'filter': True},
                                    {'field': 'seq', 'headerName': 'SEKANSLAR', 'filter': True},
                                    {'field': 'seq_len', 'headerName': 'Sekans uzunlukları', 'filter': True},
                                    {'field': 'gc', 'headerName': '%GC', 'filter': True},
                                ],
                                columnSize="sizeToFit",
                                defaultColDef={
                                    "resizable": True,
                                    "sortable": True,
                                    "filter": True,
                                    'editable': True,
                                    "minWidth": 125
                                },
                                dashGridOptions={
                                    'pagination': True,
                                    "rowSelection": "multiple",
                                    "undoRedoCellEditing": True,
                                    "undoRedoCellEditingLimit": 20,
                                    "editType": "fullRow",
                                },
                            ),

                            html.Hr(),

                            dcc.Graph(
                                figure=px.pie(
                                    data_frame=df_TABLE.to_dict("records"),
                                    names=[j for seq in df_TABLE['seq'] for j in seq],
                                    labels={'seq': "Nükleotit"},
                                ).update_traces(
                                    textposition='auto',
                                    textinfo='value+percent+label',
                                    legendgrouptitle={'text': 'Nükleotitler'}
                                )
                            ),

                            html.Hr(className="border border-danger"),
                        ],

                        style={"margin": 30},
                    )

                elif file_format == "genbank":

                    qualifiers = []

                    table = []

                    for record in records:
                        table.append({'İD': record.id, 'Tanım': record.description, 'Sekans': str(record.seq),
                                      'Sekans Uzunluğu': len(str(record.seq)), '%GC': gc_fraction(str(record.seq))})
                        if record.features:
                            for feature in record.features:
                                qualifiers.append(feature.qualifiers)

                    df_qualifiers = pd.DataFrame(i for i in qualifiers)

                    file_reading_dash_app.layout = html.Div(
                        children=[

                            html.H4("Genbank dosyası okuması", className="fw-bolder text-primary m-2 "),

                            html.Hr(className="border border-danger"),

                            dag.AgGrid(
                                id="df-table",
                                style={'width': '100%'},
                                rowData=pd.DataFrame(table).to_dict("records"),
                                columnDefs=[{"field": str(i)} for i in pd.DataFrame(table).columns],
                                columnSize="sizeToFit",
                                defaultColDef={
                                    "resizable": True,
                                    "sortable": True,
                                    "filter": True,
                                    'editable': True,
                                    "minWidth": 125
                                },
                                dashGridOptions={
                                    'pagination': True,
                                    "rowSelection": "multiple",
                                    "undoRedoCellEditing": True,
                                    "undoRedoCellEditingLimit": 20,
                                    "editType": "fullRow",
                                },
                            ),

                            html.Hr(),

                            dag.AgGrid(
                                id="df-table",
                                style={'width': '100%'},
                                rowData=df_qualifiers.to_dict("records"),
                                columnDefs=[{"field": str(i)} for i in df_qualifiers.columns],
                                columnSize="sizeToFit",
                                defaultColDef={
                                    "resizable": True,
                                    "sortable": True,
                                    "filter": True,
                                    'editable': True,
                                    "minWidth": 125
                                },
                                dashGridOptions={
                                    'pagination': True,
                                    "rowSelection": "multiple",
                                    "undoRedoCellEditing": True,
                                    "undoRedoCellEditingLimit": 20,
                                    "editType": "fullRow",
                                },
                            ),

                            html.Hr(),

                            dcc.Graph(
                                figure=px.pie(
                                    data_frame=pd.DataFrame(table).to_dict("records"),
                                    names=[j for seq in pd.DataFrame(table)['Sekans'] for j in seq],

                                ).update_traces(
                                    textposition='auto',
                                    textinfo='value+percent+label',
                                    legendgrouptitle={'text': 'Nükleotitler'},
                                )
                            ),

                            html.Hr(className="border border-danger"),

                        ], style={"margin": 30},
                    )

            finally:

                handle.close()

                file_obj.delete()

            return HttpResponseRedirect(f"/laboratory/bioinformatic/app/{file_format}-dosya-sonuc")

        else:

            form = FileReadingForm()

    return render(request, "bioinformatic/form.html", {'form': form, 'title': 'DOSYA OKUMASI'})

def fastq_stats(request):
    external_stylesheets = ['https://cdn.jsdelivr.net/npm/bootstrap@5.3.2/dist/css/bootstrap.min.css']
    external_scripts = ['https://cdn.jsdelivr.net/npm/bootstrap@5.3.2/dist/js/bootstrap.bundle.min.js']

    fastq_dash_app = DjangoDash('fastq-istatistik', external_stylesheets=external_stylesheets,
                                external_scripts=external_scripts)

    obj = FileModel.objects.get(records__user=request.user, records__tool="DOSYA OKUMA")

    handle = gzip.open(obj.file.path, 'rt', encoding='utf-8')

    records = SeqIO.parse(handle, obj.records.reading_file_format)

    cnt = defaultdict(int)

    for rec in records:
        for letter in rec.seq:
            cnt[letter] += 1

    tot = sum(cnt.values())

    df = pd.DataFrame(
        {
            'nuc': cnt.keys(),
            'count': cnt.values(),
            'per': [100 * i / tot for i in cnt.values()],
        }
    )

    columnDefs = [
        {"field": "nuc", "headerName": "Nükleotit", 'filter': True},
        {"field": "count", "headerName": "Nükleotit sayısı", 'filter': True},
        {"field": "per", "headerName": "Nükleotit Yüzdesi", 'filter': True},
    ]

    grid = dag.AgGrid(
        id="fastq-table",
        columnDefs=columnDefs,
        rowData=df.to_dict("records"),
        columnSize="sizeToFit",
        defaultColDef={
            "resizable": True,
            "sortable": True,
            "filter": True,
            'editable': True,
            "minWidth": 125
        },
        dashGridOptions={
            'pagination': True,
            "rowSelection": "multiple",
            "undoRedoCellEditing": True,
            "undoRedoCellEditingLimit": 20,
            "editType": "fullRow",
        },
    )

    fastq_dash_app.layout = html.Div(
        children=[
            html.P("FASTQ DOSYASI BİLGİLERi", className="fw-bolder text-primary mt-3"),
            grid

        ], className="m-5"
    )

    handle.close()

    obj.file.delete()

    return HttpResponseRedirect("/laboratory/bioinformatic/app/fastq-istatistik/")


def stats_view(request, user):
    global fig_lr
    if request.user.is_anonymous:
        from django.conf import settings
        messages.error(request, "Lütfen Giriş Yapınız")
        return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))

    external_stylesheets = ['https://cdn.jsdelivr.net/npm/bootstrap@5.3.2/dist/css/bootstrap.min.css']
    external_scripts = ['https://cdn.jsdelivr.net/npm/bootstrap@5.3.2/dist/js/bootstrap.bundle.min.js']

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

        df_count = pd.DataFrame({
            'Canlı': [dna.name for dna in object_list],
            'ADENİN': [a.sequence.count("A") for a in object_list],
            'TİMİN': [a.sequence.count("T") for a in object_list],
            'GUANİN': [a.sequence.count("G") for a in object_list],
            'SİTOZİN': [a.sequence.count("C") for a in object_list],
            'TOPLAM': [a.sequence.count("A") + a.sequence.count("T") + a.sequence.count("G") + a.sequence.count("C") for
                       a in object_list],
        })

        nuc_count = DjangoDash(f"nucleotit_table", external_stylesheets=external_stylesheets,
                               external_scripts=external_scripts, add_bootstrap_links=True)

        nuc_count.layout = html.Div([
            dag.AgGrid(
                style={'width': '100%'},
                rowData=df_count.to_dict("records"),
                columnSize="sizeToFit",
                defaultColDef={
                    "resizable": True,
                    "sortable": True,
                    "filter": True,
                    'editable': True,
                    "minWidth": 125
                },
                dashGridOptions={
                    'pagination': True,
                    "rowSelection": "multiple",
                    "undoRedoCellEditing": True,
                    "undoRedoCellEditingLimit": 20,
                    "editType": "fullRow",
                },
                columnDefs=[
                    {'field': 'Canlı', 'headerName': 'Canlı', 'filter': True},
                    {'field': 'ADENİN', 'headerName': 'ADENİN', 'filter': True},
                    {'field': 'TİMİN', 'headerName': 'TİMİN', 'filter': True},
                    {'field': 'GUANİN', 'headerName': 'GUANİN', 'filter': True},
                    {'field': 'SİTOZİN', 'headerName': 'SİTOZİN', 'filter': True},
                    {'field': 'TOPLAM', 'headerName': 'TOPLAM', 'filter': True},
                ]
            ),
        ])

        scatter_map = DjangoDash(f"l{request.user}r", external_stylesheets=external_stylesheets,
                                 external_scripts=external_scripts, add_bootstrap_links=True)

        scatter_map.layout = html.Div(

            [

                html.P('GRAFİK', className='fw-bolder text-primary mt-3'),

                dcc.Graph(
                    figure=px.scatter(
                        data_frame=df, x="DNA", y="PROTEİN",
                        title="DNA - PROTEİN Lineer Regrasyon",
                        trendline="ols",
                        labels={
                            "dna_seq_len": "DNA SEKANS UZUNLUĞU",
                            "pro_seq_len": "PROTEİN SEKANS UZUNLUĞU",
                            'name': 'Canlı'
                        },
                    )
                )

            ], style={'margin': 50}
        )

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

            'url': HttpResponseRedirect(f"/laboratory/bioinformatic/app/l{request.user}r/").url

        }

    else:

        return render(request, "exception/page-404.html", {'msg': "Verilere Ulaşılamadı"})

    return render(request, template_name="bioinformatic/stats/stats.html", context=context)


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
