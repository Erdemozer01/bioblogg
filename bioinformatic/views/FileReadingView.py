import os
import gzip
import subprocess, sys

import Bio.Application
from Bio import bgzf
from django.db.models import Q
from django.shortcuts import *
from django.views import generic
from django.contrib import messages
from bioinformatic.forms import FileReadingForm, TranslateForm, MultipleSeqAlignmentFileForm
from bioinformatic.models import BioinformaticModel, RecordModel, FileModel
from Bio import SeqIO
from Bio.SeqUtils import GC, gc_fraction
import pandas as pd
import plotly.express as px
from django_plotly_dash import DjangoDash
from dash import Dash, html, dcc, dash_table, Output, Input
from pathlib import Path
from Bio.Seq import Seq
import dash_ag_grid as dag
from collections import defaultdict
from Bio.Align.Applications import ClustalwCommandline, ClustalOmegaCommandline, MuscleCommandline
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceCalculator
from Bio import Phylo
import dash_cytoscape as cyto
from bioinformatic.generate_tree import generate_elements
import dash_bio as dashbio

BASE_DIR = Path(__file__).resolve().parent.parent.parent


def PhylogeneticTree(request):
    global clustal_omega_exe, muscle_cline, tree_alg, elements, stylesheet, layout, styles, tree, clustalw2_exe, stats, stats_file
    stats = ""
    if request.user.is_anonymous:
        from django.conf import settings
        messages.error(request, "Lütfen Giriş Yapınız")
        return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))

    external_stylesheets = ['https://cdn.jsdelivr.net/npm/bootstrap@5.3.2/dist/css/bootstrap.min.css']
    external_scripts = ['https://cdn.jsdelivr.net/npm/bootstrap@5.3.2/dist/js/bootstrap.bundle.min.js']

    if BioinformaticModel.objects.filter(user=request.user, tool="Filogenetik Ağaç").exists():
        BioinformaticModel.objects.filter(user=request.user, tool="Filogenetik Ağaç").delete()

    form = MultipleSeqAlignmentFileForm(request.POST or None, request.FILES or None)

    if request.method == "POST":

        if request.FILES['file'].size > 25 * 1024 * 1024:
            messages.error(request, "Dosya boyutu 25mb dan fazladır.")
            return HttpResponseRedirect(request.path)

        if form.is_valid():
            tree_obj = BioinformaticModel.objects.create(user=request.user, tool="Filogenetik Ağaç")
            in_file = tree_obj.records_files.create(file=request.FILES['file'])
            out_file = tree_obj.records_files.create(file=f"{request.user}_alignment.fasta")
            xml_tree_file = tree_obj.records_files.create(file=f"{request.user}_tree.xml")
            cluster_file = tree_obj.records_files.create(file=f"{request.user}_cluster.csv")
            scores_file = tree_obj.records_files.create(file=f"{request.user}_scores.txt")

            tools = form.cleaned_data['alignment_tools']

            tree_alg = form.cleaned_data['tree_alg']

            tree_app = DjangoDash(f'{tree_alg}',
                                  external_stylesheets=external_stylesheets,
                                  external_scripts=external_scripts,
                                  add_bootstrap_links=True,
                                  title=f'{tree_alg} AĞACI OLUŞTRMA'.upper()
                                  )
            try:

                if tools == "muscle":

                    if sys.platform.startswith('win32'):
                        muscle_cline = os.path.join(BASE_DIR, 'bioinformatic', 'programs',
                                                    "muscle3.8.31_i86win32.exe")
                    elif sys.platform.startswith('linux'):
                        muscle_cline = os.path.join(BASE_DIR, 'bioinformatic', 'programs',
                                                    'muscle3.8.425_i86linux32')

                    muscle_cline_tool = MuscleCommandline(
                        muscle_cline,
                        input=in_file.file.path,
                        out=out_file.file.path,

                    )

                    if sys.platform.startswith('win32'):
                        assert os.path.isfile(muscle_cline)
                        stdout, stderr = muscle_cline_tool()

                    elif sys.platform.startswith('linux'):
                        muscle_exe = os.path.join(BASE_DIR, 'bioinformatic', 'programs', 'muscle3.8.425_i86linux32')
                        muscle_result = subprocess.check_output(
                            [muscle_exe, "-in", in_file.file.path, "-out", out_file.file.path])

                elif tools == "clustalw2":
                    stats_file = tree_obj.records_files.create(file=f"{request.user}_stats.txt")
                    if sys.platform.startswith('win32'):
                        clustalw2_exe = os.path.join(BASE_DIR, 'bioinformatic', 'programs', 'clustalw2.exe')

                    elif sys.platform.startswith('linux'):
                        clustalw2_exe = os.path.join(BASE_DIR, 'bioinformatic', 'programs', 'clustalw2')

                    clustalw_cline = ClustalwCommandline(
                        clustalw2_exe,
                        infile=in_file.file.path,
                        outfile=out_file.file.path,
                        align=True,
                        outorder="ALIGNED",
                        convert=True,
                        output="FASTA",
                        stats=stats_file.file.path,
                        clustering=tree_alg,
                        score="percent"
                    )

                    assert os.path.isfile(clustalw2_exe)
                    stdout, stderr = clustalw_cline()

            except Bio.Application.ApplicationError:
                messages.error(request, "Uygulama Hatası! Hatalı dosya türü girdiniz.")
                tree_obj.delete()
                return HttpResponseRedirect(request.path)

            alignment = AlignIO.read(handle=out_file.file.path, format='fasta')

            data = open(out_file.file.path, 'r', encoding="utf-8").read()

            calculator = DistanceCalculator('identity')

            constructor = DistanceTreeConstructor(calculator, method=tree_alg)

            trees = constructor.build_tree(alignment)

            Phylo.write(trees=trees, file=xml_tree_file.file.path, format="phyloxml")

            tree = Phylo.read(xml_tree_file.file.path, 'phyloxml')

            try:

                nodes, edges = generate_elements(tree)

                elements = nodes + edges

            except ZeroDivisionError:
                messages.error(request, "Hata! Protein dosyası girmiş olabilirsiniz")
                return HttpResponseRedirect(request.path)

            layout = {"name": "preset"}

            stylesheet = [
                {
                    "selector": ".nonterminal",
                    "style": {
                        "label": "data(confidence)",
                        "background-opacity": 0,
                        "text-halign": "left",
                        "text-valign": "top",
                    },
                },
                {"selector": ".support", "style": {"background-opacity": 0}},
                {
                    "selector": "edge",
                    "style": {
                        "source-endpoint": "inside-to-node",
                        "target-endpoint": "inside-to-node",
                    },
                },
                {
                    "selector": ".terminal",
                    "style": {
                        "label": "data(name)",
                        "width": 10,
                        "height": 10,
                        "text-valign": "center",
                        "text-halign": "right",
                        "background-color": "#222222",
                    },
                },
            ]

            tree_app.layout = html.Div(
                [

                    html.P("FİLOGENETİK AĞAÇ OLUŞTURMA", className="text-primary fw-bolder mb-4"),

                    html.P(f"{tools} aracıyla {tree_alg} ağacı oluşturuldu".upper(),
                           className="text-primary fw-bolder mb-4"),

                    cyto.Cytoscape(
                        id='cytoscape-usage-phylogeny',
                        elements=elements,
                        stylesheet=stylesheet,
                        layout=layout,
                        style={
                            'height': '95vh',
                            'width': '100%'
                        }
                    ),

                    html.Hr(),

                    html.P("Alignment Haritası", className="fw-bolder"),

                    dashbio.AlignmentChart(
                        id='my-default-alignment-viewer',
                        data=data,
                        height=900,
                        tilewidth=30,
                    ),

                ], className="m-5"
            )

            @tree_app.callback(
                Output('cytoscape-usage-phylogeny', 'stylesheet'),
                Input('cytoscape-usage-phylogeny', 'mouseoverEdgeData')
            )
            def color_children(edgeData):
                if edgeData is None:
                    return stylesheet

                if 's' in edgeData['source']:
                    val = edgeData['source'].split('s')[0]
                else:
                    val = edgeData['source']

                children_style = [{
                    'selector': 'edge[source *= "{}"]'.format(val),
                    'style': {
                        'line-color': 'blue'
                    }
                }]

                return stylesheet + children_style

        return HttpResponseRedirect(f"/laboratory/bioinformatic/app/{tree_alg}")

    return render(request, 'bioinformatic/form.html', {'form': form, 'title': 'Filogenetik ağaç oluşturma'})


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

        if request.FILES['file'].size > 25 * 1024 * 1024:
            messages.error(request, "Dosya boyutu en fazla 25mb olmalıdır.")
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

            records = SeqIO.parse(handle, file_format)

            file_reading_dash_app = DjangoDash(
                name=f"{file_format}-dosya-sonuc",
                external_stylesheets=external_stylesheets,
                external_scripts=external_scripts,
                title=f"{file_format} dosyası okuması sonuçları".upper()
            )

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

            elif file_format == "fastq":
                n_cnt = defaultdict(int)

                cnt_qual = defaultdict(int)

                qual_pos = defaultdict(list)

                handle = gzip.open(file_obj.file.path, 'rt', encoding='utf-8')

                try:
                    records = SeqIO.parse(handle, file_format)

                    for rec in records:
                        for i, letter in enumerate(rec.seq):
                            pos = i + 1
                            if letter == 'N':
                                n_cnt[pos] += 1
                        for pos, quality in enumerate(rec.letter_annotations['phred_quality']):
                            cnt_qual[quality] += 1
                            if pos < 25 or quality == 40:
                                continue
                            pos = pos + 1
                            qual_pos[pos].append(quality)

                except gzip.BadGzipFile:

                    handle = open(file_obj.file.path, 'r', encoding='utf-8')
                    records = SeqIO.parse(handle, file_format)

                    for rec in records:
                        for i, letter in enumerate(rec.seq):
                            pos = i + 1
                            if letter == 'N':
                                n_cnt[pos] += 1
                        for pos, quality in enumerate(rec.letter_annotations['phred_quality']):
                            cnt_qual[quality] += 1
                            if pos < 25 or quality == 40:
                                continue
                            pos = pos + 1
                            qual_pos[pos].append(quality)

                tot = sum(cnt_qual.values())
                seq_len = max(n_cnt.keys())
                positions = range(1, seq_len + 1)

                df_table = pd.DataFrame([[qual, 100. * cnt / tot, cnt] for qual, cnt in cnt_qual.items()],
                                        columns=['KALİTE SKORLARI', 'YÜZDE', 'SAYISI'])

                file_reading_dash_app.layout = html.Div(
                    [

                        html.H4(f'{file_format} DOSYASI SONUÇLARI'.upper(), className="text-primary"),

                        html.Hr(),

                        html.Button("Tabloyu indir", id="csv-button", n_clicks=0,
                                    className="btn btn-outline-primary btn-sm float-sm-end mb-1"),

                        html.P("TABLO"),

                        dag.AgGrid(
                            style={'width': '100%'},
                            id="export-data-grid",
                            rowData=df_table.to_dict("records"),
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
                            columnDefs=[{'field': f'{col}', 'headerName': f'{col}', 'filter': True} for col in
                                        df_table.columns],
                            csvExportParams={
                                "fileName": "table_data.csv",
                            },
                        ),

                        html.Hr(),

                        html.P("Grafik Seçiniz", className="fw-bolder ml-1"),

                        dcc.Dropdown(
                            id="select",
                            options={
                                "distribution_n": "N dağılım grafiği",
                                "phred_score": "Phred skor dağılım grafiği",
                            },
                        ),

                        html.Div(id="figure-output", className="mt-3"),

                    ], className="m-5"
                )

                @file_reading_dash_app.callback(
                    Output("figure-output", "children"),
                    Input("select", "value"),
                )
                def figure(select):
                    if select == "distribution_n":
                        return dcc.Graph(figure=px.line(x=positions, y=[n_cnt[x] for x in positions]))
                    elif select == "phred_score":
                        vps = []
                        poses = list(qual_pos.keys())
                        poses.sort()
                        for pos in poses:
                            vps.append(qual_pos[pos])
                        return dcc.Graph(
                            figure=px.box(data_frame=vps, y=poses,
                                          labels={'variable': 'Nükleotit Pozisyonu', 'value': 'PHERED skorları'},
                                          title='PHERED skorları dağılımı')
                        )

                @file_reading_dash_app.callback(
                    Output("export-data-grid", "exportDataAsCsv"),
                    Input("csv-button", "n_clicks"),
                    prevent_initial_call=True,
                )
                def export_data_as_csv(n_clicks):
                    if n_clicks:
                        return True
                    return False

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
