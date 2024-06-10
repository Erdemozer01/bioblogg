import os
import gzip
import subprocess, sys
import Bio.Application
import parmed
from django.shortcuts import *
from django.contrib import messages
from bioinformatic.forms import FileReadingForm, MultipleSeqAlignmentFileForm, BlastForm, MoleculeViewForm
from bioinformatic.models import BioinformaticModel
from Bio import SeqIO, SearchIO
from Bio.SeqUtils import gc_fraction
import pandas as pd
import plotly.express as px
from django_plotly_dash import DjangoDash
from dash import Dash, html, dcc, dash_table, Output, Input, State
import dash_bootstrap_components as dbc
from pathlib import Path
from Bio.Blast import NCBIWWW, NCBIXML
import dash_ag_grid as dag
from collections import defaultdict
from Bio.Align.Applications import ClustalwCommandline, ClustalOmegaCommandline, MuscleCommandline
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceCalculator
from Bio import Phylo
import dash_cytoscape as cyto
from bioinformatic.generate_tree import generate_elements
import dash_bio as dashbio
from dash_bio.utils import PdbParser, create_mol3d_style
import dash_bio.utils.ngl_parser as ngl_parser

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
                        clustering=tree_alg,
                        score="percent"
                    )

                    assert os.path.isfile(clustalw2_exe)
                    stdout, stderr = clustalw_cline()

                elif tools == "omega":
                    if sys.platform.startswith('win32'):
                        clustal_omega_exe = os.path.join(BASE_DIR, 'bioinformatic', 'programs',
                                                         'clustal-omega-1.2.2-win64/clustalo.exe')
                    elif sys.platform.startswith('linux'):
                        clustal_omega_exe = os.path.join(BASE_DIR, 'bioinformatic', 'programs',
                                                         'clustalo-1.2.4-Ubuntu-32-bit')

                    clustal_omega_cline = ClustalOmegaCommandline(
                        clustal_omega_exe,
                        infile=in_file.file.path,
                        outfile=out_file.file.path,
                        force=True,
                        verbose=True,
                        auto=True,
                        usekimura="yes"
                    )

                    if sys.platform.startswith('win32'):
                        assert os.path.isfile(clustal_omega_exe)
                        stdout, stderr = clustal_omega_cline()

                    elif sys.platform.startswith('linux'):
                        subprocess.Popen(str(clustal_omega_cline), stdin=subprocess.PIPE,
                                         stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE, universal_newlines=True,
                                         shell=(sys.platform != "win32"))

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

                    dbc.NavbarSimple(
                        children=[
                            html.A(
                                [dbc.NavItem(dbc.NavLink("Biyoinformatik", className="text-white"))],
                                href=redirect("bioinformatic:home").url
                            ),

                            html.A(
                                [dbc.NavItem(dbc.NavLink("Geri", className="text-white")), ],
                                href=redirect("bioinformatic:pyhlo_tree").url,
                            ),

                        ],
                        brand="Filogenetik ağaç oluşturma".title(),
                        brand_href=str(request.path),
                        color="primary",
                        sticky="top",
                    ),

                    html.P(f"{tools} aracıyla {tree_alg} ağacı oluşturuldu".upper(),
                           className="text-primary fw-bolder mb-4 mt-3"),

                    cyto.Cytoscape(
                        id='cytoscape-usage-phylogeny',
                        elements=elements,
                        stylesheet=stylesheet,
                        layout=layout,
                        style={
                            'height': '80vh',
                            'width': '90%'
                        }, className="mx-auto"
                    ),

                    html.Hr(),

                    html.P("Alignment Haritası", className="fw-bolder"),

                    dashbio.AlignmentChart(
                        id='my-default-alignment-viewer',
                        data=data,
                        height=1000,
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

    return render(request, 'bioinformatic/form.html', {'form': form, 'title': 'Filogeneni ve Alignment Haritası'})


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


def blast(request):
    global result_handle
    if request.user.is_anonymous:
        from django.conf import settings
        messages.error(request, "Lütfen Giriş Yapınız")
        return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))

    external_stylesheets = ['https://cdn.jsdelivr.net/npm/bootstrap@5.3.2/dist/css/bootstrap.min.css']
    external_scripts = ['https://cdn.jsdelivr.net/npm/bootstrap@5.3.2/dist/js/bootstrap.bundle.min.js']

    app = DjangoDash('blast', external_stylesheets=external_stylesheets,
                     external_scripts=external_scripts, add_bootstrap_links=True,
                     title='BLAST')

    if BioinformaticModel.objects.filter(user=request.user, tool="BLAST").exists():
        BioinformaticModel.objects.filter(user=request.user, tool="BLAST").delete()

    form = BlastForm(request.POST or None, request.FILES or None)

    blast_obj = BioinformaticModel.objects.create(
        user=request.user,
        tool="BLAST",
    )

    if request.method == "POST":
        if form.is_valid():

            type = form.cleaned_data['type']
            program = form.cleaned_data['program']
            database = form.cleaned_data['database']

            if type == 'file':
                file = blast_obj.records_files.create(file=request.FILES['file'])
                record = next(SeqIO.parse(file.file.path, format="fasta"))
                result_handle = NCBIWWW.qblast(program=program, database=database, sequence=record.seq)
            elif type == 'gi':
                gi = form.cleaned_data['gi']
                result_handle = NCBIWWW.qblast(str(program), str(database), str(gi))

            blast_xml = blast_obj.records_files.create(file=f"{request.user}_blast.xml")
            result_file = blast_obj.records_files.create(file="result.txt")

            with open(blast_xml.file.path, "w") as out_handle:
                out_handle.write(result_handle.read())

            blast_qresult = SearchIO.read(blast_xml.file.path, "blast-xml")

            blast_records = NCBIXML.parse(open(blast_xml.file.path, 'r'))

            file_obj = open(result_file.file.path, "w")

            for hsp in blast_qresult:
                for blast_record in blast_records:
                    for alignment in blast_record.alignments:
                        for hsps in alignment.hsps:
                            file_obj.writelines(f"{hsp}" + 3 * "\n")
                            file_obj.writelines(f"{hsps}" + 3 * "\n")

            app.layout = dbc.Container(
                [
                    html.H1("Blast Sonuçları", className="mt-3"),
                    html.Hr(),
                    dbc.Row(
                        [

                            html.Button("Sonuçları indir", id="btn-download-blast",
                                        style={'float': 'right', 'marginBottom': '20px'},
                                        className='btn btn-primary col-4'),
                            dcc.Download(id="download-blast"),

                            dcc.Textarea(id='align-textarea',
                                         value=open(result_file.file.path, 'r').read(),
                                         readOnly=True,
                                         style={'width': '100%', 'height': '300px'}),
                        ],
                        align="center", className="m-4"
                    ),
                ],

            )

            @app.callback(
                Output("download-blast", "data"),
                Input("btn-download-blast", "n_clicks"),
                Input('align-textarea', 'value'),
                prevent_initial_call=True,
            )
            def func(n_clicks, values):
                return dict(content=values, filename=f"{request.user}_blast.txt")

            file_obj.close()

            blast_obj.delete()

            return HttpResponseRedirect("/laboratory/bioinformatic/app/blast/")

    return render(request, 'bioinformatic/form.html', {'form': form, 'title': 'BLAST'})


def molecule_viewer(request):
    if request.user.is_anonymous:
        from django.conf import settings
        messages.error(request, "Lütfen Giriş Yapınız")
        return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))

    if BioinformaticModel.objects.filter(user=request.user, tool="molecule_view").exists():
        BioinformaticModel.objects.filter(user=request.user, tool="molecule_view").delete()

    external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

    app = DjangoDash('3d_molecule_view', external_stylesheets=external_stylesheets,

                     title='3D MOLEKÜL GÖRÜNTÜLEME')

    form = MoleculeViewForm(request.POST or None, request.FILES or None)

    if request.method == "POST":

        if form.is_valid():

            if request.FILES['file'].size > 25 * 1024 * 1024:
                messages.error(request, "Dosya boyutu en fazla 25mb olmalıdır.")
                return HttpResponseRedirect(request.path)

            file = form.cleaned_data["file"]

            obj = BioinformaticModel.objects.create(
                user=request.user,
                tool="molecule_view",
            )

            file_obj = obj.records_files.create(file=file)

            try:

                parser = PdbParser(file_obj.file.path)

                data = parser.mol3d_data()

                styles = create_mol3d_style(
                    data['atoms'], visualization_type='cartoon', color_element='residue'
                )

                df = pd.DataFrame(data["atoms"])
                df = df.drop_duplicates(subset=['residue_name'])
                df['positions'] = df['positions'].apply(lambda x: ', '.join(map(str, x)))

            except ValueError:
                messages.error(request, "Sayfayı yenilediğiniz İçin veriler kaybolmuştur")

                return redirect('bioinformatic:3d_molecule_view')

            except parmed.exceptions.FormatNotFound:

                messages.error(request, "Hatalı Dosya Formatı")

                return redirect('bioinformatic:3d_molecule_view')

            columns = [
                {'name': 'Seri', 'id': 'serial'},
                {'name': 'Adı', 'id': 'name'},
                {'name': 'ELEMENT', 'id': 'elem'},
                {'name': 'POZİSYON', 'id': 'positions'},
                {'name': 'Kütle Büyüklüğü', 'id': 'mass_magnitude'},
                {'name': 'İNDEX', 'id': 'residue_index'},
                {'name': 'Bölge Adı', 'id': 'residue_name'},
                {'name': 'Zincir', 'id': 'chain'}
            ]

            app.layout = html.Div(
                [

                    html.H4("3D MOLEKÜL GÖRÜNTÜLEME", className="mt-3"),

                    html.P(f'İD : {file.name[:4]},'),

                    html.Hr(),

                    dash_table.DataTable(
                        id="zooming-specific-residue-table",
                        columns=columns,
                        data=df.to_dict("records"),
                        row_selectable="single",
                        page_size=10,
                        filter_action='native',
                        filter_options={"placeholder_text": "Filtrele"},
                        editable=True,
                        style_cell={'textAlign': 'center'},
                        style_header={
                            'backgroundColor': 'white',
                            'fontWeight': 'bold'
                        },
                    ),
                    html.Hr(),
                    dashbio.Molecule3dViewer(
                        id="zooming-specific-molecule3d-zoomto",
                        modelData=data,
                        styles=styles,
                        style={'marginRight': 'auto', 'marginLeft': 'auto'},
                        width=950, height=600
                    ),
                ], className="container"
            )

            @app.callback(
                Output("zooming-specific-molecule3d-zoomto", "zoomTo"),
                Output("zooming-specific-molecule3d-zoomto", "labels"),
                Input("zooming-specific-residue-table", "selected_rows"),
                prevent_initial_call=True
            )
            def residue(selected_row):
                row = df.iloc[selected_row]
                row['positions'] = row['positions'].apply(lambda x: [float(x) for x in x.split(',')])

                return [
                    {
                        "sel": {"chain": row["chain"], "resi": row["residue_index"]},
                        "animationDuration": 1500,
                        "fixedPath": True,
                    },

                    [
                        {
                            "text": "Bölge ADI: {}".format(row["residue_name"].values[0]),
                            "position": {
                                "x": row["positions"].values[0][0],
                                "y": row["positions"].values[0][1],
                                "z": row["positions"].values[0][2],
                            },
                        }
                    ],

                ]

        return HttpResponseRedirect("/laboratory/bioinformatic/app/3d_molecule_view/")

    return render(request, 'bioinformatic/form.html', {'form': form, 'title': '3D MOLECULE GÖRÜNTÜLEME'})
