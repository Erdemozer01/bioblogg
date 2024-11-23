##### MEHMET ERDEM ÖZER, mozer232@posta.pau.edu.tr ######
from pathlib import Path
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import xGC_skew
import pandas as pd
import plotly.express as px
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from dash import html, dcc, Input, Output, State, Patch
from django.shortcuts import *
from django_plotly_dash import DjangoDash
import dash_ag_grid as dag
import plotly.figure_factory as ff
from Bio import Align
from Bio.Align import substitution_matrices
import dash_bootstrap_components as dbc
import dash_daq as daq
from Bio.Emboss.Applications import Primer3Commandline
from bioinformatic.primer_design import design_primers

BASE_DIR = Path(__file__).resolve().parent.parent.parent


def alignment_score(request):
    external_stylesheets = [dbc.themes.BOOTSTRAP]
    alignment_score = DjangoDash('alignment_score', external_stylesheets=external_stylesheets,
                                 title="ALİGNMENT SKOR", add_bootstrap_links=True)

    alignment_score.layout = html.Div(

        [
            dbc.NavbarSimple(
                children=[
                    dbc.NavItem(dbc.NavLink("Blog", href=HttpResponseRedirect(
                        reverse("blog:anasayfa")).url, external_link=True)),
                    dbc.DropdownMenu(
                        children=[
                            dbc.DropdownMenuItem("Biyoinformatik",
                                                 href=HttpResponseRedirect(reverse("bioinformatic:home")).url,
                                                 external_link=True),
                            dbc.DropdownMenuItem("Biyoistatislik",
                                                 href=HttpResponseRedirect(reverse("biostatistic:home")).url,
                                                 external_link=True),
                            dbc.DropdownMenuItem("Coğrafi Bilgi sistemleri",
                                                 href=HttpResponseRedirect(reverse("cbs")).url,
                                                 external_link=True),
                            dbc.DropdownMenuItem("Laboratuvarlar",
                                                 href=HttpResponseRedirect(reverse("lab_home")).url,
                                                 external_link=True),
                        ],
                        nav=True,
                        in_navbar=True,
                        label="Laboratuvarlar",

                    ),
                ],
                brand="ALİGNMENT SCORE",
                brand_href=HttpResponseRedirect(reverse("bioinformatic:alignment_score")).url,
                color="primary",
                dark=True,
                brand_external_link=True,
                sticky='top',
                className="shadow-lg p-3 bg-body rounded container"
            ),

            html.Div(
                [
                    dbc.Container(

                        [

                            dbc.Row(
                                [
                                    dbc.Col(
                                        [
                                            html.P("Alignment Mod", className='text-primary fw-bolder'),

                                            dcc.Dropdown(
                                                id="mod",
                                                options=[
                                                    {'label': 'LOCAL', 'value': 'local'},
                                                    {'label': 'GLOBAL', 'value': 'global'},

                                                ], value='local',
                                            ),

                                            html.P("Alignment MATRİX", style={'marginTop': '10px'},
                                                   className='text-primary fw-bolder'),

                                            dcc.Dropdown(id="matrix", options={
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

                                            html.P("Hedef sekans", style={'marginTop': '10px'},
                                                   className='text-primary fw-bolder'),

                                            dcc.Textarea(id="target_seq", placeholder="Sekans Giriniz",
                                                         style={'marginRight': '10px', 'width': '100%',
                                                                'height': '100px'},
                                                         className="form-control"),

                                            html.P("Sorgu sekans", style={'marginTop': '10px'},
                                                   className='text-primary fw-bolder'),

                                            dcc.Textarea(id="query_seq", placeholder="Sekans Giriniz",
                                                         style={'marginRight': '10px', 'width': '100%',
                                                                'height': '100px'},
                                                         className='form-control'),
                                        ], md=4
                                    ),
                                    dbc.Col(
                                        [
                                            html.Div(id="output"),
                                        ], md=8
                                    )
                                ])

                        ], className="shadow-lg p-3 mb-5 bg-body rounded mt-2", fluid=False
                    ),
                ],
            ),
        ], style={'marginTop': "2%"}
    )

    @alignment_score.callback(
        Output("output", "children"),
        Input("mod", "value"),
        Input("matrix", "value"),
        Input("target_seq", "value"),
        Input("query_seq", "value"),
    )
    def update_output(mod, matrix, target_seq, query_seq):
        target_seq = str(target_seq).replace(" ", '')
        target_seq = str(target_seq).replace("\n", "")
        target_seq = str(target_seq).replace("\t", "")
        target_seq = str(target_seq).replace("\r", "")
        target_seq = str(target_seq).replace("0", "")
        target_seq = str(target_seq).replace("1", "")
        target_seq = str(target_seq).replace("2", "")
        target_seq = str(target_seq).replace("3", "")
        target_seq = str(target_seq).replace("4", "")
        target_seq = str(target_seq).replace("5", "")
        target_seq = str(target_seq).replace("6", "")
        target_seq = str(target_seq).replace("7", "")
        target_seq = str(target_seq).replace("8", "")
        target_seq = str(target_seq).replace("9", "")
        target_seq = str(target_seq).upper()

        query_seq = str(query_seq).replace(" ", "")
        query_seq = str(query_seq).replace("\n", "")
        query_seq = str(query_seq).replace("\t", "")
        query_seq = str(query_seq).replace("\r", "")
        query_seq = str(query_seq).replace("0", "")
        query_seq = str(query_seq).replace("1", "")
        query_seq = str(query_seq).replace("2", "")
        query_seq = str(query_seq).replace("3", "")
        query_seq = str(query_seq).replace("4", "")
        query_seq = str(query_seq).replace("5", "")
        query_seq = str(query_seq).replace("6", "")
        query_seq = str(query_seq).replace("7", "")
        query_seq = str(query_seq).replace("8", "")
        query_seq = str(query_seq).replace("9", "")
        query_seq = str(query_seq).upper()

        aligner = Align.PairwiseAligner()

        aligner.mode = str(mod)

        aligner.substitution_matrix = substitution_matrices.load(str(matrix))

        alignments = aligner.align(target_seq, query_seq)

        alignment = alignments[0]

        count = alignment.counts()

        values = f'Boşluk : {count.gaps}  Benzerlik : {count.identities},  Eşleşmeyen: {count.mismatches}\n\nMatriks: {matrix}\n\n {alignment.substitutions}\nMOD: {str(mod).upper()}\n\nScore: {alignments.score}\n\nAlignment:\n\n' + str(
            alignment)

        return html.Div(
            [

                html.P("SONUÇLAR", className='fw-bolder mt-1 text-primary'),

                dcc.Textarea(
                    id='align-textarea',
                    value=f"Hedef sekans uzunluğu: {len(target_seq)}, "
                          f"%GC: {gc_fraction(target_seq) * 100}, "
                          f"Sorgu sekans uzunluğu: {len(query_seq)}, "
                          f"%GC: {gc_fraction(query_seq) * 100}\n\n" + values,

                    readOnly=True,
                    style={'width': '100%', 'height': '400px'},
                    disabled=True
                ),

                html.Button("Sonuçları indir", id="btn-download-align",
                            className='btn btn-primary mr-2 mt-2 col-12'),

                dcc.Download(id="download-align"),

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

    return HttpResponseRedirect("/laboratuvar/bioinformatic/app/alignment_score/")


def sequence_analiz(request):
    external_stylesheets = [dbc.themes.ZEPHYR]
    external_scripts = ["https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/js/bootstrap.bundle.min.js"]

    sequence_analiz = DjangoDash("sequence_analiz", external_stylesheets=external_stylesheets,
                                 title='Sekans Analiz', add_bootstrap_links=True, external_scripts=external_scripts)

    sequence_analiz.layout = html.Div(

        [

            dbc.NavbarSimple(
                children=[
                    dbc.NavItem(dbc.NavLink("Blog", href=HttpResponseRedirect(
                        reverse("blog:anasayfa")).url, external_link=True)),
                    dbc.DropdownMenu(
                        children=[
                            dbc.DropdownMenuItem("Biyoinformatik",
                                                 href=HttpResponseRedirect(reverse("bioinformatic:home")).url,
                                                 external_link=True),
                            dbc.DropdownMenuItem("Biyoistatislik",
                                                 href=HttpResponseRedirect(reverse("biostatistic:home")).url,
                                                 external_link=True),
                            dbc.DropdownMenuItem("Coğrafi Bilgi sistemleri",
                                                 href=HttpResponseRedirect(reverse("cbs")).url,
                                                 external_link=True),
                            dbc.DropdownMenuItem("Laboratuvarlar",
                                                 href=HttpResponseRedirect(reverse("lab_home")).url,
                                                 external_link=True),
                        ],
                        nav=True,
                        in_navbar=True,
                        label="Laboratuvarlar",

                    ),
                ],
                brand="SEKANS ANALİZİ",
                brand_href=HttpResponseRedirect(reverse("bioinformatic:dna_seq_read")).url,
                color="primary",
                dark=True,
                brand_external_link=True,
                sticky='top',
                className="shadow-lg p-3 bg-body rounded container"
            ),

            html.Div(
                dbc.Container(
                    [
                        dbc.Row([
                            dbc.Col(
                                [

                                    html.Label("Nükleotit tipi", className="fw-bolder mt-2"),
                                    dcc.Dropdown(
                                        id='nuc_type',
                                        options=[
                                            {'label': 'DNA', 'value': 'dna'},
                                            {'label': 'RNA', 'value': 'rna'},
                                        ],
                                        value='dna'
                                    ),

                                    html.Label("Grafik tipi", className="fw-bolder mt-2"),
                                    dcc.Dropdown(
                                        id='grap_type',
                                        options=[
                                            {'label': 'Dağılım', 'value': 'dist'},
                                            {'label': 'Histogram', 'value': 'hist'},
                                        ],
                                        value='dist'
                                    ),

                                    html.Label("Sekans", className="fw-bolder mt-2"),
                                    dcc.Textarea(id="sequence", placeholder="Sekans Giriniz",
                                                 className="form-control mb-1", style={'height': 200}),
                                ], md=4
                            ),

                            dbc.Col(
                                [

                                    html.Div(id="seq_output"),
                                    html.Hr(),
                                    html.Div(id="fig_output", className="shadow-lg rounded mb-1"),
                                ], md=8
                            ),

                        ]),

                    ], className="shadow-lg p-4 mb-4 bg-body rounded mt-2", fluid=False
                ),
            ),

        ], className="mt-4",
    )

    @sequence_analiz.callback(
        Output("seq_output", "children"),
        Output("fig_output", "children"),
        Input("sequence", "value"),
        Input("nuc_type", "value"),
        Input("grap_type", "value"),
    )
    def update_output(sequence, nuc_type, grap_type):
        if sequence:
            if nuc_type == "dna":
                sequence = sequence.upper()
                ig = []

                for i in sequence:
                    if not i in ['A', 'C', 'G', 'T']:
                        ig.append(i)

                for i in ig:
                    if i in sequence:
                        sequence = sequence.replace(i, "")

            elif nuc_type == "rna":
                sequence = sequence.upper()
                ig = []
                for i in sequence:
                    if not i in ['A', 'C', 'G', 'U']:
                        ig.append(i)

                for i in ig:
                    if i in sequence:
                        sequence = sequence.replace(i, "")

                if not "U" in sequence:
                    return html.P(["DNA SEKANSI OLABİLİR..."]), ""

            nuc_df = pd.DataFrame({
                'positions': [pos for pos in range(len(sequence))],
                'nucleotit': [seq for seq in sequence]
            })

            a = nuc_df[(nuc_df['nucleotit']) == "A"]
            t = nuc_df[(nuc_df['nucleotit']) == "T"]
            g = nuc_df[(nuc_df['nucleotit']) == "G"]
            c = nuc_df[(nuc_df['nucleotit']) == "C"]

            seq_results = (
                f"SEKANS UZUNLUĞU: {len(sequence)}, %GC: {gc_fraction(sequence)}, A: {sequence.count('A')}, "
                f"T: {sequence.count('T')}, G: {sequence.count('G')}, C: {sequence.count('C')}")

            dist_figure = dcc.Graph(figure=ff.create_distplot(
                [a['positions'], t['positions'], g['positions'], c['positions']],
                ['A', 'T', 'G', 'C'],
                bin_size=[.1, .25, .5, 1], show_hist=False,
                curve_type="normal",
            ).update_layout({'title': 'Nükleotit Dağılım Grafiği'}))

            hist_grap = dcc.Graph(figure=px.histogram(
                nuc_df.to_dict("records"), x="nucleotit", color="nucleotit",
                title="Nükleotit Sayıları",
            ).update_yaxes(title="Nükletotit Sayısı").update_xaxes(title="Nükleotit"))

            if grap_type == "dist":
                return seq_results, dist_figure

            elif grap_type == "hist":

                return seq_results, hist_grap

        else:
            return html.P("Henüz bir sekans dizisi girmediniz ! ".upper(),
                          className='text-danger mx-auto text-center mt-3')

    return HttpResponseRedirect("/laboratuvar/bioinformatic/app/sequence_analiz/")


def translation(request):
    external_stylesheets = [dbc.themes.BOOTSTRAP]
    translate_app = DjangoDash("translate_dna", external_stylesheets=external_stylesheets,
                               add_bootstrap_links=True, title='PROTEİN SENTEZİ')

    translate_app.layout = dbc.Card(

        [

            dbc.NavbarSimple(
                children=[
                    dbc.NavItem(dbc.NavLink("Blog", href=HttpResponseRedirect(
                        reverse("blog:anasayfa")).url, external_link=True)),

                    dbc.DropdownMenu(
                        children=[
                            dbc.DropdownMenuItem("Biyoinformatik",
                                                 href=HttpResponseRedirect(reverse("bioinformatic:home")).url,
                                                 external_link=True),
                            dbc.DropdownMenuItem("Biyoistatislik",
                                                 href=HttpResponseRedirect(reverse("biostatistic:home")).url,
                                                 external_link=True),
                            dbc.DropdownMenuItem("Coğrafi Bilgi sistemleri",
                                                 href=HttpResponseRedirect(reverse("cbs")).url,
                                                 external_link=True),
                            dbc.DropdownMenuItem("Laboratuvarlar",
                                                 href=HttpResponseRedirect(reverse("lab_home")).url,
                                                 external_link=True),
                        ],
                        nav=True,
                        in_navbar=True,
                        label="Laboratuvarlar"

                    ),
                ],
                brand="Protein Sentezi",
                brand_href=HttpResponseRedirect(reverse("bioinformatic:dna_seq_translate")).url,
                color="primary",
                dark=True,
                brand_external_link=True,
                sticky='top',
                className="shadow-lg p-3 bg-body rounded"
            ),

            dbc.CardBody([
                dbc.Row([
                    dbc.Col(

                        [

                            dbc.Tabs([
                                dbc.Tab(
                                    label="Sekans",

                                    children=[

                                        dcc.Textarea(id="seq", placeholder="Sekans Giriniz",
                                                     className="form-control mb-1 mt-3", style={'height': 200}),

                                        html.Label("Sekans başlangıcı", className="fw-bolder mt-2"),

                                        dbc.Input(id="start", type="number", placeholder="Sekans başlangıcı", min=0,
                                                  className="mr-1"),

                                        html.Label("Sekans bitişi", className="fw-bolder mt-2"),

                                        dbc.Input(id="stop", type="number", placeholder="Sekans bitişi", min=0,
                                                  className="mr-1"),

                                        html.Label("İstenmeyen sekans dizisi", className="fw-bolder mt-2"),

                                        dbc.Input(id="discard", type="text", placeholder="İstenmeyen sekans dizisi",
                                                  className="mr-1"),
                                    ]
                                ),

                                dbc.Tab(
                                    label="Dönüşüm tablosu",

                                    children=[
                                        html.Label("Dönüşüm Tablosu", className="fw-bolder mt-2"),

                                        dcc.Dropdown(id="table", options={
                                            "1": "Standart Kod",
                                            "2": "Omurgalı Mitokondri Kodu",
                                            "3": "Maya Mitokondri Kodu",
                                            "4": "Küf, Protozoon ve Kölenterat Mitokondri Kodu ve Mikoplazma/Spiroplasma Kodu",
                                            "5": "Omurgasız Mitokondri Kodu",
                                            "6": "Siliat, Dasikladas ve Heksamita Nükleer Kodu",
                                            "9": "Ekinoderm ve Yassı Solucan Mitokondri Kodu",
                                            "10": "Euplotid Nükleer Kodu",
                                            "11": "Bakteriyel, Arkeal ve Bitki Plastid Kodu, prokaryotik virüsler",
                                            "12": "Alternatif Maya Nükleer Kodu",
                                            "13": "Ascidian Mitokondri Kodu",
                                            "14": "Alternatif Yassı Solucan Mitokondri Kodu",
                                            "16": "Klorofis Mitokondri Kodu",
                                            "21": "Trematod Mitokondriyal Kodu",
                                            "22": "Scenedesmus obliquus Mitokondri Kodu",
                                            "23": "Thraustochytrium Mitokondri Kodu",
                                            "24": "Rhabdopleuridae Mitokondri Kodu",
                                            "25": "Aday Bölüm SR1 ve Gracilibacteria Kodu",
                                            "26": "Pachysolen tannophilus Nükleer Kodu",
                                            "27": "Karyorelict Nükleer Kodu",
                                            "28": "Kondilostoma Nükleer Kodu",
                                            "29": "Mezodinyum Nükleer Kodu",
                                            "30": "Peritrich Nükleer Kodu",
                                            "31": "Blastocrithidia Nükleer Kodu",
                                            "33": "Cephalodiscidae Mitokondriyal UAA-Tyr Kodu"
                                        }, value="1", searchable=True),

                                        html.Small(" ***Dönüşüm tablosu verileri ", className='txt-bold'),

                                        html.A(' NCBI ',
                                               href="https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi",
                                               target="_blank",
                                               style={'text-decoration': 'none'}),

                                        html.Small(' sitesiden alınmıştır.'),

                                        html.Div([
                                            daq.BooleanSwitch(
                                                id="to_stop",
                                                on=False,
                                                label="Stop Kodonları",
                                                labelPosition="top",
                                                className="fw-bolder float-left mt-2",
                                            )
                                        ]),

                                    ]
                                )
                            ])

                        ], md=4, className="mt-3"
                    ),

                    dbc.Col(
                        [

                            dbc.Tabs([
                                dbc.Tab(
                                    label="Sonuçlar",

                                    children=[

                                        html.P(id="seq_inf", className="mt-2"),

                                        html.P(id="protein_inf"),

                                        html.Label("Komplement", style={'font-weight': 'bold'}),

                                        dcc.Textarea(
                                            id="complement",
                                            className="form-control",
                                            readOnly=True
                                        ),

                                        html.Label("Reverse Komplement", style={'font-weight': 'bold'}),

                                        dcc.Textarea(
                                            id="reverse_complement",
                                            className="form-control",
                                            readOnly=True
                                        ),

                                        html.Label("Transcribe", style={'font-weight': 'bold'}),

                                        dcc.Textarea(
                                            id="transcribe",
                                            className="form-control",
                                            readOnly=True
                                        ),

                                        html.Label(f"Protein Sekansı", style={'font-weight': 'bold'}),

                                        dcc.Textarea(
                                            id="protein",
                                            className="form-control",
                                            readOnly=True
                                        ),
                                    ]
                                )
                            ]),
                        ], md=8, className="mt-3"
                    )
                ]),
            ]),

        ], className="shadow-lg p-3 mb-5 bg-body rounded mt-2 mr-2 ml-2"
    )

    @translate_app.callback(
        Output("seq_inf", "children"),
        Output("protein_inf", "children"),
        Output("complement", "value"),
        Output("reverse_complement", "value"),
        Output("transcribe", "value"),
        Output("protein", "value"),

        Input("seq", "value"),
        Input("start", "value"),
        Input("stop", "value"),
        Input("discard", "value"),
        Input("table", "value"),
        Input("to_stop", "on"),
    )
    def update_seq(seq, start, stop, discard, table, to_stop):

        sequence = str(seq).upper()
        ig = []

        for i in sequence:
            if not i in ['A', 'C', 'G', 'T']:
                ig.append(i)

        for i in ig:
            if i in sequence:
                sequence = sequence.replace(i, "")

        if discard:
            sequence = sequence.replace(str(discard).upper(), "")

        if start or stop:
            sequence = sequence[start:stop]

        if to_stop is True:
            trans_seq = Seq(sequence).translate(table=int(table))
        else:
            trans_seq = Seq(sequence).translate(table=int(table)).replace("*", "")

        seq_inf = (f"Girilen sekans uzunluğu: {len(sequence)}, %GC: {gc_fraction(sequence)}, A: {sequence.count('A')} "
                   f"T: {sequence.count('T')}, G: {sequence.count('G')}, C: {sequence.count('C')} \n"
                   )
        protein_inf = f"Protein Sekans Uzunluğu: {len(trans_seq)}, Stop kodonları sayısı: {trans_seq.count('*')}"

        complement = str(Seq(sequence).complement())

        reverse_complement = str(Seq(sequence).reverse_complement())

        transcribe = str(Seq(sequence).transcribe())

        return seq_inf, protein_inf, complement, reverse_complement, transcribe, str(trans_seq)

    return HttpResponseRedirect("/laboratuvar/bioinformatic/app/translate_dna/")


def Kmer_SeqSlicing(request):
    external_stylesheets = [dbc.themes.BOOTSTRAP]
    external_scripts = ["https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/js/bootstrap.bundle.min.js"]

    seq_input = DjangoDash('kmer_seq_slicing', external_stylesheets=external_stylesheets,
                           title='Kmer oluşturma', add_bootstrap_links=True, external_scripts=external_scripts)

    seq_input.layout = html.Div(

        [

            dbc.NavbarSimple(
                children=[
                    dbc.NavItem(dbc.NavLink("Blog", href=HttpResponseRedirect(
                        reverse("blog:anasayfa")).url, external_link=True)),
                    dbc.DropdownMenu(
                        children=[
                            dbc.DropdownMenuItem("Biyoinformatik",
                                                 href=HttpResponseRedirect(
                                                     reverse("bioinformatic:home")).url,
                                                 external_link=True),
                            dbc.DropdownMenuItem("Biyoistatislik",
                                                 href=HttpResponseRedirect(
                                                     reverse("biostatistic:home")).url,
                                                 external_link=True),
                            dbc.DropdownMenuItem("Coğrafi Bilgi sistemleri",
                                                 href=HttpResponseRedirect(reverse("cbs")).url,
                                                 external_link=True),
                            dbc.DropdownMenuItem("Laboratuvarlar",
                                                 href=HttpResponseRedirect(
                                                     reverse("lab_home")).url,
                                                 external_link=True),
                        ],
                        nav=True,
                        in_navbar=True,
                        label="Laboratuvarlar",
                        className="float-right",

                    ),
                ],
                brand="KMER OLUŞTURMA",
                brand_href=HttpResponseRedirect(reverse("bioinformatic:dna_seq_slice")).url,
                color="primary",
                dark=True,
                brand_external_link=True,
                sticky='top',
                className="shadow-lg bg-body rounded mt-1 mb-1 ",
            ),

            dbc.Card(
                [
                    dbc.CardBody(
                        dbc.Row(
                            [
                                dbc.Col(
                                    [

                                        html.P("DNA Sekans giriniz", className='fw-bolder'),

                                        html.Div(
                                            dcc.Textarea(
                                                id="seq",
                                                placeholder="Sekans Giriniz",
                                                style={'width': '100%', 'height': '200px'},
                                                className="form-control"
                                            ), style={'whiteSpace': 'pre-line'}
                                        ),

                                        dbc.Label(["Sekans başlangıcı"], className="fw-bolder mt-1"),
                                        dcc.Input(id="start", type="number", placeholder="Sekans başlangıcı", min=0,
                                                  value=0,
                                                  className="form-control mb-1"),

                                        dbc.Label(["Sekans bitişi"], className="fw-bolder"),
                                        dcc.Input(id="stop", type="number", placeholder="Sekans bitişi", min=1,
                                                  value=1, className="form-control mb-1"),

                                        dbc.Label(["Kmer Uzunluğu"], className="fw-bolder"),
                                        dcc.Input(id="kmer", type="number", placeholder="Kmer Uzunluğu", min=1,
                                                  className="form-control mb-1"),

                                        dbc.Label(["İstenmeyen sekans dizisi"], className="fw-bolder"),
                                        dcc.Input(id="discard", type="text", placeholder="İstenmeyen sekans dizisi",
                                                  className="form-control mb-1"),
                                    ], md=4
                                ),
                                dbc.Col(
                                    [
                                        html.Div(id="output"),
                                    ], md=8
                                ),
                            ],

                        )
                    )
                ]

            ),

        ], className="shadow-lg p-3 bg-body rounded mr-1 ml-1 mt-1"
    )

    @seq_input.callback(
        Output("output", "children"),
        Input("seq", "value"),
        Input("start", "value"),
        Input("stop", "value"),
        Input("kmer", "value"),
        Input("discard", "value"),
    )
    def update_output(seq, start, stop, kmer, discard):
        sequence = seq.upper()
        ig = []

        for i in sequence:
            if not i in ['A', 'C', 'G', 'T']:
                ig.append(i)

        for i in ig:
            if i in sequence:
                sequence = sequence.replace(i, "")

        if discard:
            sequence = sequence.replace(str(discard).upper(), "")

        if start or stop:
            sequence = sequence[start:stop]

        kmer_list = []

        def getKmers(sequence, size):
            for x in range(0, len(sequence) - size):
                yield sequence[x:x + size]

        for km in getKmers(sequence, kmer):
            kmer_list.append(km)

        df_kmer = pd.DataFrame(
            {
                'id': [i + 1 for i in range(len(kmer_list))],
                'kmer': kmer_list,
                '%gc': [gc_fraction(gc) for gc in kmer_list],
            }
        )

        return html.Div(
            [

                html.H4("SONUÇLAR", className='fw-bolder mt-2'),

                html.P(f"SEKANS UZUNLUĞU: {len(sequence)}, %GC: {gc_fraction(sequence)}",
                       style={'marginTop': '10px'}, className="mt-2"),

                html.P("İstenilen Sekans", className='fw-bolder'),

                dcc.Textarea(value=sequence, style={'width': '100%', 'height': '100px'}, className="form-control"),

                html.P("Kmerler", className="fw-bolder mt-1"),

                dag.AgGrid(
                    style={'width': '100%'},
                    rowData=df_kmer.to_dict("records"),
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
                        {'field': 'id', 'headerName': 'İD', 'filter': True},
                        {'field': 'kmer', 'headerName': 'KMERS', 'filter': True},
                        {'field': '%gc', 'headerName': '%GC', 'filter': True},
                    ]
                ),
            ],
        ), html.Br(), html.Hr(),

    return HttpResponseRedirect("/laboratuvar/bioinformatic/app/kmer_seq_slicing/")


def create_frame_seq(request):
    external_stylesheets = ['https://cdn.jsdelivr.net/npm/bootstrap@5.3.2/dist/css/bootstrap.min.css']
    external_scripts = ['https://cdn.jsdelivr.net/npm/bootstrap@5.3.2/dist/js/bootstrap.bundle.min.js']

    frame_seq_app = DjangoDash('create_frame_seq', external_stylesheets=external_stylesheets,
                               external_scripts=external_scripts, add_bootstrap_links=True,
                               title='Çerçeve sekans oluşturma'.upper())

    frame_seq_app.layout = html.Div(

        [

            ## NAVBAR ##
            dbc.NavbarSimple(
                children=[
                    dbc.NavItem(dbc.NavLink("Blog", href=HttpResponseRedirect(
                        reverse("blog:anasayfa")).url, external_link=True)),
                    dbc.DropdownMenu(
                        children=[
                            dbc.DropdownMenuItem("Biyoinformatik",
                                                 href=HttpResponseRedirect(
                                                     reverse("bioinformatic:home")).url,
                                                 external_link=True),
                            dbc.DropdownMenuItem("Biyoistatislik",
                                                 href=HttpResponseRedirect(
                                                     reverse("biostatistic:home")).url,
                                                 external_link=True),
                            dbc.DropdownMenuItem("Coğrafi Bilgi sistemleri",
                                                 href=HttpResponseRedirect(reverse("cbs")).url,
                                                 external_link=True),
                            dbc.DropdownMenuItem("Laboratuvarlar",
                                                 href=HttpResponseRedirect(
                                                     reverse("lab_home")).url,
                                                 external_link=True),
                        ],
                        nav=True,
                        in_navbar=True,
                        label="Laboratuvarlar",
                        className="float-right",

                    ),
                ],
                brand="Çerçeve Sekans Oluşturma",
                brand_href=HttpResponseRedirect(reverse("bioinformatic:create_frame_seq")).url,
                color="primary",
                dark=True,
                brand_external_link=True,
                sticky='top',
                className="shadow-lg bg-body rounded mt-2 mb-1 ",
            ),

            dbc.Card(
                dbc.CardBody(
                    dbc.Row(
                        [
                            dbc.Col(
                                [

                                    dbc.Label("Sekans Giriniz", className="fw-bolder mt-1"),
                                    dcc.Textarea(
                                        id="seq", placeholder="Sekans Giriniz",
                                        className="form-control mb-2", style={'height': 100}
                                    ),

                                    dbc.Label("Nükleotit pozisyonu", className="fw-bolder mt-1"),
                                    dcc.Input(id="nuc_pos", type="number", placeholder="Nükleotit pozisyonu", min=0,
                                              className="form-control mb-2", value=1, ),

                                    dbc.Label("Çerçeve sekans uzunluğu", className="fw-bolder mt-1"),
                                    dcc.Input(id="frame_seq_len", type="number", placeholder="Çerçeve sekans uzunluğu",
                                              min=10, className="form-control mb-2", value=50, ),

                                    dbc.Label("Sekans başlangıcı", className="fw-bolder mt-1"),
                                    dcc.Input(id="start", type="number", placeholder="Sekans başlangıcı", min=0,
                                              className="form-control mb-2", value=0, ),

                                    dbc.Label("Sekans bitişi", className="fw-bolder mt-1"),
                                    dcc.Input(id="stop", type="number", placeholder="Sekans bitişi", min=2,
                                              className="form-control mb-1"),

                                    dbc.Label("İstenmeyen sekans dizisi", className="fw-bolder mt-2"),
                                    dcc.Input(id="discard", type="text", placeholder="İstenmeyen sekans dizisi",
                                              className="form-control mb-2"),
                                ], md=4
                            ),

                            dbc.Col(
                                [
                                    html.Div(id="output_seq"),
                                ], md=8
                            ),
                        ]
                    )
                )
            ),

        ], className="shadow-lg p-3 bg-body rounded mr-1 ml-1 mt-1"
    )

    @frame_seq_app.callback(
        Output("output_seq", "children"),

        Input("seq", "value"),
        Input("nuc_pos", "value"),
        Input("frame_seq_len", "value"),
        Input("start", "value"),
        Input("stop", "value"),
        Input("discard", "value"),
    )
    def update_output(seq, nuc_pos, frame_seq_len, start, stop, discard):

        sequence = str(seq).upper()

        ig = []

        for i in sequence:
            if not i in ['A', 'C', 'G', 'T']:
                ig.append(i)

        for i in ig:
            if i in sequence:
                sequence = sequence.replace(i, "")

        if discard:
            sequence = sequence.replace(str(discard).upper(), "")

        print(sequence)

        if start or stop:
            sequence = sequence[start:stop]

        df = pd.DataFrame({
            'seq': [k for k in sequence]
        }, index=[i + 1 for i in range(len(sequence))])

        nuc_position = df.to_dict()['seq'].get(nuc_pos)

        df_frame = pd.DataFrame(
            {
                'seq_len': [len(sequence[:50] + str(nuc_position) + sequence[50:frame_len]) for frame_len in
                            range(50, frame_seq_len)],
                'seq_frame': [sequence[:50] + ' ' + str(nuc_position) + ' ' + sequence[50:frame_len] for frame_len in
                              range(50, frame_seq_len)]
            }
        )

        return html.Div([
            html.P(
                f"SEKANS UZUNLUĞU: {len(sequence)}, %GC: {gc_fraction(sequence)}, Nükleotid Pozisyonu : {nuc_position}",
                style={'marginTop': '20px'}

            ),

            html.P("Çerçeve Sekanslar", className="fw-bolder"),

            html.Div([
                dag.AgGrid(
                    id="frame_seq_table",
                    style={'width': '100%'},
                    rowData=df_frame.to_dict("records"),
                    columnSize="sizeToFit",
                    defaultColDef={"resizable": True, "sortable": True, "filter": True, 'editable': True,
                                   "minWidth": 125},
                    dashGridOptions={'pagination': True, "rowSelection": "multiple"},
                    columnDefs=[
                        {'field': 'seq_len', 'headerName': 'UZUNLUK', 'filter': True},
                        {'field': 'seq', 'headerName': 'SEKANS', 'filter': True},
                    ]
                ),
                html.Div(id="output_frame"),
            ]),

            html.Hr(),
        ]),

    @frame_seq_app.callback(
        Output("output_frame", "children"),
        Input("frame_seq_table", "cellClicked"),
        prevent_initial_call=True,
    )
    def display_cell_clicked_on(cell):
        if cell is None:
            return html.P("Çerçeve sekans seçimi henüz yapılmadı!", className="text-danger")

        return html.Div(
            [
                dcc.Textarea(
                    value=f"Çerçeve Sekans: \n{cell['value']}",
                    style={'width': '100%', 'height': '100px'}),
            ]
        ),

    return HttpResponseRedirect("/laboratuvar/bioinformatic/app/create_frame_seq/")


def TemperatureMeltingView(request):
    external_stylesheets = [dbc.themes.BOOTSTRAP]

    mt_app = DjangoDash('temp-melt', external_stylesheets=external_stylesheets,
                        title='Primer Erime Sıcaklığı', add_bootstrap_links=True)

    mt_app.layout = html.Div(

        [
            dbc.NavbarSimple(
                children=[
                    dbc.NavItem(dbc.NavLink("Blog", href=HttpResponseRedirect(
                        reverse("blog:anasayfa")).url, external_link=True)),
                    dbc.DropdownMenu(
                        children=[
                            dbc.DropdownMenuItem("Biyoinformatik",
                                                 href=HttpResponseRedirect(reverse("bioinformatic:home")).url,
                                                 external_link=True),
                            dbc.DropdownMenuItem("Biyoistatislik",
                                                 href=HttpResponseRedirect(reverse("biostatistic:home")).url,
                                                 external_link=True),
                            dbc.DropdownMenuItem("Coğrafi Bilgi sistemleri",
                                                 href=HttpResponseRedirect(reverse("cbs")).url,
                                                 external_link=True),
                            dbc.DropdownMenuItem("Laboratuvarlar",
                                                 href=HttpResponseRedirect(reverse("lab_home")).url,
                                                 external_link=True),
                        ],
                        nav=True,
                        in_navbar=True,
                        label="Laboratuvarlar",

                    ),
                ],
                brand="Primer Sıcaklık Hesaplama",
                brand_href=HttpResponseRedirect(reverse("bioinformatic:temp_melt")).url,
                color="primary",
                dark=True,
                brand_external_link=True,
                sticky='top',
                className="shadow-lg p-3 bg-body rounded-2 container"
            ),

            html.Div(
                [
                    dbc.Container(

                        [

                            dbc.Row(
                                [
                                    dbc.Col(
                                        [
                                            dcc.Tabs(
                                                id='mol3d-tabs', children=[

                                                    dcc.Tab(
                                                        label='Sekans',
                                                        children=html.Div(
                                                            className='control-tab mt-2',
                                                            children=[

                                                                html.Label("Metodlar", className="fw-bolder mt-1"),
                                                                dcc.Dropdown(
                                                                    id='method-type',
                                                                    placeholder="Seçiniz",
                                                                    options=[
                                                                        {'label': 'Klasik', 'value': 'Tm_Wallace'},
                                                                        {'label': 'GC içerik tabanlı',
                                                                         'value': 'Tm_GC'},
                                                                        {'label': 'Termodinamik Tabanlı (Önerilen)',
                                                                         'value': 'Tm_NN'},
                                                                    ], value="Tm_NN"
                                                                ),

                                                                dbc.Label("Sekans", className="fw-bolder mt-1"),
                                                                dcc.Textarea(placeholder="Sekans", id="seq",
                                                                             style={'height': 100},
                                                                             className="form-control"),

                                                                dbc.Label("Komplement Sekans",
                                                                          className="fw-bolder mt-1"),
                                                                dcc.Textarea(placeholder="Komplement Sekans",
                                                                             id="c_seq",
                                                                             style={'height': 100},
                                                                             className="form-control"),
                                                            ]
                                                        )
                                                    ),

                                                    dcc.Tab(
                                                        label='Konsantrasyon',
                                                        children=[
                                                            html.Div(
                                                                [

                                                                    html.P(
                                                                        [
                                                                            "Aşagıdaki değerler von Ahsen et al., 2001 referans alınarak ayarlanmıştır."],
                                                                    ),

                                                                    dbc.Label("Sodyum", className="fw-bolder mt-1"),
                                                                    dbc.Input(placeholder="Sodyum", type="number",
                                                                              id="Na", value=50, min=0,
                                                                              className="form-control"),

                                                                    dbc.Label("Potayum", className="fw-bolder mt-1"),
                                                                    dbc.Input(placeholder="Potasyum", type="number",
                                                                              id="K", value=0, min=0,
                                                                              className="form-control"),

                                                                    dbc.Label("Tris", className="fw-bolder mt-1"),
                                                                    dbc.Input(placeholder="Tris", type="number",
                                                                              id="Tris", value=0, min=0,
                                                                              className="form-control"),

                                                                    dbc.Label("Magnezyum", className="fw-bolder mt-1"),
                                                                    dbc.Input(placeholder="Magnezyum", type="number",
                                                                              id="Mg", value=0, min=0,
                                                                              className="form-control"),

                                                                    dbc.Label("dNTPs", className="fw-bolder mt-1"),
                                                                    dbc.Input(placeholder="dNTPs", type="number",
                                                                              id="dNTPs", value=0, min=0,
                                                                              className="form-control"),

                                                                    dbc.Label("Tuz", className="fw-bolder mt-1"),
                                                                    dbc.Input(placeholder="Tuz", type="number",
                                                                              id="saltcorr", value=5, min=0,
                                                                              max=7, className="form-control"),

                                                                ],
                                                                className="mx-auto"
                                                            ),
                                                        ]
                                                    ),

                                                ],
                                            ),

                                        ], md=8,
                                    ),

                                    dbc.Col(
                                        [
                                            html.Div(id="output"),
                                        ], md=4
                                    )
                                ])

                        ], className="shadow-lg p-4 bg-body rounded-2 mt-1", fluid=False
                    ),
                ],
            ),

        ], style={'marginTop': "1%"}
    )

    @mt_app.callback(
        Output("output", 'children'),
        # seq
        Input("method-type", "value"),
        Input("seq", "value"),
        Input("c_seq", "value"),
        # kons
        Input("Na", "value"),
        Input("K", "value"),
        Input("Tris", "value"),
        Input("Mg", "value"),
        Input("dNTPs", "value"),
        Input("saltcorr", "value"),
    )
    def temp_melting(method, seq, c_seq, Na, K, Tris, Mg, dNTPs, saltcorr):
        global temp_melting

        if seq is not None:

            primer_len = len(seq)

            if method == "Tm_Wallace":
                temp_melting = mt.Tm_Wallace(seq=seq)

            elif method == "Tm_GC":
                temp_melting = mt.Tm_GC(seq=seq, Na=Na, K=K, Tris=Tris, Mg=Mg, dNTPs=dNTPs)

            elif method == "Tm_NN":

                temp_melting = mt.Tm_NN(seq=seq, c_seq=c_seq, Na=Na, K=K, Tris=Tris, Mg=Mg, dNTPs=dNTPs,
                                        saltcorr=saltcorr)

            return html.Div(
                [
                    html.P(f"Sekans Uzunluğu : {primer_len}"),
                    html.P(f"Erime Sıcaklığı : {temp_melting}"),
                ], className="text-primary"
            )

        else:
            return html.P(["Sekans girilmedi"], className="text-danger")

    return HttpResponseRedirect("/laboratuvar/bioinformatic/app/temp-melt/")


def PrimerDesign(request):
    external_stylesheets = [dbc.themes.BOOTSTRAP]
    app = DjangoDash('primer-design', external_stylesheets=external_stylesheets,
                     title="Primer Dizaynı", add_bootstrap_links=True)

    app.layout = html.Div(
        [
            dbc.NavbarSimple(
                children=[
                    dbc.NavItem(dbc.NavLink("Blog", href=HttpResponseRedirect(
                        reverse("blog:anasayfa")).url, external_link=True)),
                    dbc.DropdownMenu(
                        children=[
                            dbc.DropdownMenuItem("Biyoinformatik",
                                                 href=HttpResponseRedirect(reverse("bioinformatic:home")).url,
                                                 external_link=True),
                            dbc.DropdownMenuItem("Biyoistatislik",
                                                 href=HttpResponseRedirect(reverse("biostatistic:home")).url,
                                                 external_link=True),
                            dbc.DropdownMenuItem("Coğrafi Bilgi sistemleri",
                                                 href=HttpResponseRedirect(reverse("cbs")).url,
                                                 external_link=True),
                            dbc.DropdownMenuItem("Laboratuvarlar",
                                                 href=HttpResponseRedirect(reverse("lab_home")).url,
                                                 external_link=True),
                        ],
                        nav=True,
                        in_navbar=True,
                        label="Laboratuvarlar",

                    ),
                ],
                brand="Primer Dizaynı",
                brand_href=HttpResponseRedirect(reverse("bioinformatic:primer_design")).url,
                color="primary",
                dark=True,
                brand_external_link=True,
                sticky='top',
                className="shadow-lg p-3 bg-body rounded-2 container"
            ),

            html.Div(
                [
                    dbc.Container(

                        [

                            dbc.Row(
                                [
                                    dbc.Col(
                                        [
                                            dbc.Tabs(
                                                id='primer-design',
                                                children=[
                                                    dcc.Tab(
                                                        label='Sekans',
                                                        children=html.Div(
                                                            className='control-tab mt-2',
                                                            children=[

                                                                dbc.Label("Sekans Türü", className="fw-bolder mt-1"),

                                                                dcc.Dropdown(
                                                                    id="mod",
                                                                    options=[
                                                                        {'label': 'DNA', 'value': 'DNA'},
                                                                        {'label': 'RNA', 'value': 'RNA'},
                                                                    ], value='DNA',
                                                                ),

                                                                dbc.Label("Sekans", className="fw-bolder mt-1"),
                                                                dcc.Textarea(id="seq",
                                                                             style={'height': 200},
                                                                             className="form-control"),

                                                                html.Div(id="seq_len"),

                                                                html.Small(
                                                                    "*Aradığınız yada ekleme çıkarma yapmak istediğiniz bölgeyi CTRL + F ile tarayıcıdan arayıp sekans alanından silebilirsiniz. Sonuçlarınız otomatik olarak değişecektir.",
                                                                    className="text-danger fst-italic ")

                                                            ]
                                                        )
                                                    ),



                                                    dcc.Tab(
                                                        label='Primer',
                                                        children=[
                                                            html.Div(
                                                                [

                                                                    dbc.Label("Primer Uzunluğu",
                                                                              className="fw-bolder mt-1"),
                                                                    dbc.Input(type="number",
                                                                              id="primer_length", value=20, min=0,
                                                                              className="form-control"),

                                                                    dbc.Label("Min %GC",
                                                                              className="fw-bolder mt-1"),
                                                                    dbc.Input(type="number",
                                                                              id="gc_min", min=0, value=40,
                                                                              className="form-control"),

                                                                    dbc.Label("Max %GC", className="fw-bolder mt-1"),
                                                                    dbc.Input(type="number",
                                                                              id="gc_max", value=60, min=0,
                                                                              className="form-control"),

                                                                    dbc.Label("Min Erime Sıcaklığı (°C)",
                                                                              className="fw-bolder mt-1"),
                                                                    dbc.Input(type="number",
                                                                              id="tm_min", value=50, min=0,
                                                                              className="form-control"),

                                                                    dbc.Label("Max Erime Sıcaklığı (°C)",
                                                                              className="fw-bolder mt-1"),
                                                                    dbc.Input(type="number",
                                                                              id="tm_max", value=60, min=0,
                                                                              className="form-control"),

                                                                ],
                                                                className="mx-auto"
                                                            ),
                                                        ]
                                                    ),

                                                ],


                                            ),

                                        ], md=4,
                                    ),

                                    dbc.Col(
                                        [

                                            dbc.Tabs(
                                                dbc.Tab(
                                                    label="Sonuçlar",
                                                    children=[
                                                        html.Div(id="filter_div"),
                                                        html.Div(id="results_div"),
                                                        html.Div(id="primer_results_len", className="mt-5 mb-2"),
                                                        html.Div(id="uyarı")
                                                    ]
                                                )
                                            ),

                                        ], md=8
                                    )
                                ])

                        ], className="shadow-lg p-4 bg-body rounded-2 mt-2", fluid=False
                    ),
                ],
            ),

        ], style={'marginTop': "1%"}
    )

    @app.callback(
        Output("seq", "placeholder"),
        Input("mod", "value"),
    )
    def nuc_type(type):
        return f"{type} Sekansı Giriniz"

    @app.callback(
        Output("seq_len", "children"),
        Output("filter_div", "children"),
        Output("results_div", "children"),
        Output("primer_results_len", "children"),
        Output("uyarı", "children"),

        Input("seq", "value"),
        Input("primer_length", "value"),
        Input("gc_min", "value"),
        Input("gc_max", "value"),
        Input("tm_min", "value"),
        Input("tm_max", "value"),
        Input("mod", "value"),
    )
    def primer_design(seq, primer_length, gc_min, gc_max, tm_min, tm_max, mod):
        aligner = Align.PairwiseAligner()
        if seq:
            seq = seq.upper()
            ig = []

            for i in seq:
                if not i in ['A', 'C', 'G', 'T']:
                    ig.append(i)

            for i in ig:
                if i in seq:
                    seq = seq.replace(i, "")


            seq_len = html.Small(f"Sekans Uzunluğu : {len(seq)}nt", className="fw-bold fst-italic")

            forward_primers = []
            reverse_primers = []


            for s in range(len(seq)):
                if primer_length == len(seq[s:s + primer_length]):
                    forward_primers.append(seq[s:s + primer_length])


            for r in forward_primers[::-1]:
                if mod == "DNA":
                    reverse_primers.append(Seq(r).reverse_complement())
                elif mod == "RNA":
                    reverse_primers.append(Seq(r).reverse_complement_rna())

            last_fp = []
            last_rp = []

            for fp in forward_primers:
                if gc_min <= gc_fraction(fp) * 100 <= gc_max and tm_min <= mt.Tm_Wallace(fp) <= tm_max:
                    last_fp.append(fp)

            for rp in reverse_primers:
                if gc_min <= gc_fraction(rp) * 100 <= gc_max and tm_min <= mt.Tm_Wallace(rp) <= tm_max:
                    last_rp.append(rp)

            df = pd.DataFrame(
                {
                    'İleri Primer': [str(i) for i in last_fp],
                    'İleri Primer %GC': [gc_fraction(str(i)) * 100 for i in last_fp],
                    'İleri Primer TM(°C)': [mt.Tm_Wallace(str(i)) for i in last_fp],

                    'Geri primer': [str(i) for i in last_rp],
                    'Geri Primer %GC': [gc_fraction(str(i)) * 100 for i in last_rp],
                    'Geri Primer TM(°C)': [mt.Tm_Wallace(str(i)) for i in last_rp],
                    'Hit score': [aligner.align(Seq(f), Seq(r)).score for f, r in zip(last_fp, last_rp)],
                }
            )

            filter = dbc.Input(id="quick-filter-input", placeholder="Tabloda Ara...",
                               className="form-control float-right mb-1 mt-2")

            results = dag.AgGrid(
                id="results_table",
                columnDefs=[{"field": i, 'editable': True} for i in df.columns],
                rowData=df.to_dict("records"),
                columnSize="sizeToFit",
                dashGridOptions={"animateRows": True},
            )

            created_primer = html.P(f"{len(last_fp)} primer oluşturuldu.")

            uyarı = html.Small("*Hit puanı primerlerin eşleşme scorudur.")

            return seq_len, filter, results, created_primer, uyarı

        else:

            seq_len = html.P(f"Henüz sekans girmediniz!", className="text-danger fst-italic")

            return seq_len, None, None, None, None

    @app.callback(
        Output("results_table", "dashGridOptions"),
        Input("quick-filter-input", "value"),
        prevent_initial_call=True
    )
    def update_filter(filter_value):
        newFilter = Patch()
        newFilter['quickFilterText'] = filter_value
        return newFilter



    return HttpResponseRedirect("/laboratuvar/bioinformatic/app/primer-design/")
