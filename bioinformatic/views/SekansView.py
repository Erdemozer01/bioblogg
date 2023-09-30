from pathlib import Path
from Bio.SeqUtils import MeltingTemp as mt
import pandas as pd
import plotly.express as px
from Bio.Seq import Seq
from Bio.SeqUtils import GC, gc_fraction
from dash import html, dcc, Input, Output
from django.shortcuts import *
from django_plotly_dash import DjangoDash
import dash_ag_grid as dag
import plotly.figure_factory as ff
from Bio import Align
from Bio.Align import substitution_matrices

BASE_DIR = Path(__file__).resolve().parent.parent.parent


def alignment_score(request):
    external_stylesheets = ['https://cdn.jsdelivr.net/npm/bootstrap@5.3.2/dist/css/bootstrap.min.css']
    external_scripts = ['https://cdn.jsdelivr.net/npm/bootstrap@5.3.2/dist/js/bootstrap.bundle.min.js']

    alignment_score = DjangoDash('alignment_score', external_stylesheets=external_stylesheets,
                                 external_scripts=external_scripts, add_bootstrap_links=True, title="ALİGNMENT SKOR")

    alignment_score.layout = html.Div(

        [

            html.Div(
                [

                    html.H4('ALİGNMENT SKOR', className='text-primary fw-bolder'),

                    html.A('BİYOİNFORMATİK ANASAYFA', href=HttpResponseRedirect(reverse("bioinformatic:home")).url,
                           style={'float': 'right'}),

                    html.P("Alignment Mod", className='text-success fw-bolder'),

                    dcc.Dropdown(id="mod", options={
                        'local': 'LOCAL',
                        'global': 'GLOBAL'
                    }, searchable=True),

                    html.P("Alignment MATRİX", style={'marginTop': '10px'}, className='text-success fw-bolder'),

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

                    html.P("Hedef sekans", style={'marginTop': '10px'}, className='text-success fw-bolder'),

                    dcc.Textarea(id="seq1", placeholder="Sekans Giriniz",
                                 style={'marginRight': '10px', 'width': '100%', 'height': '100px'}),

                    html.P("Sorgu sekans", style={'marginTop': '10px'}, className='text-success fw-bolder'),

                    dcc.Textarea(id="seq2", placeholder="Sekans Giriniz",
                                 style={'marginRight': '10px', 'width': '100%', 'height': '100px'}),

                ], style={'margin': 20},
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

        values = f'Boşluk : {count.gaps}  Benzerlik : {count.identities},  Eşleşmeyen: {count.mismatches},\n\nMatriks: {mat},\n\n {alignment.substitutions}\nMOD: {str(mod).upper()}\n\nScore: {alignments.score}\n\nAlignment:\n\n' + str(
            alignment)

        return html.Div(
            [

                html.P(
                    f"Hedef sekans uzunluğu: {len(seq1)}, %GC: {GC(seq1)}, Sorgu sekans uzunluğu: {len(seq2)}, %GC: {GC(seq2)}",
                    style={'marginTop': '20px'}
                ),

                html.Hr(),

                html.P("Alignments sonuçları", style={'marginBottom': '-35px'}, className='fw-bolder'),

                html.Button("Sonuçları indir", id="btn-download-align",
                            style={'float': 'right', 'marginBottom': '20px'},
                            className='btn btn-primary mr-2'),
                dcc.Download(id="download-align"),

                html.Div(
                    [
                        dcc.Textarea(id='align-textarea', value=values, readOnly=True,
                                     style={'width': '100%', 'height': '500px'}),
                    ]
                ),

                html.Hr(),
                html.Br(),

            ], id='output', style={'margin': 20}

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


def sequence_analiz(request):
    external_stylesheets = ['https://cdn.jsdelivr.net/npm/bootstrap@5.3.2/dist/css/bootstrap.min.css']
    external_scripts = ['https://cdn.jsdelivr.net/npm/bootstrap@5.3.2/dist/js/bootstrap.bundle.min.js']
    sequence_analiz = DjangoDash("sequence_analiz", external_stylesheets=external_stylesheets,
                                 external_scripts=external_scripts, add_bootstrap_links=True, title='Sekans Analiz')

    sequence_analiz.layout = html.Div(

        [

            html.H4('Sekans Analiz', className='text-primary'),

            html.A('BİYOİNFORMATİK ANASAYFA', href=HttpResponseRedirect(reverse("bioinformatic:home")).url,
                   style={'float': 'right'}),

            html.P("Sekans", className='fw-bolder'),

            html.Div(
                [
                    dcc.Textarea(id="sekans", placeholder="Sekans Giriniz",
                                 className="form-control mb-1", style={'height': 200}),
                    dcc.Input(id="start", type="number", placeholder="Sekans başlangıcı", min=0,
                              style={'marginLeft': '2px'}),
                    dcc.Input(id="stop", type="number", placeholder="Sekans bitişi", min=0,
                              style={'marginLeft': '2px'}),
                    dcc.Input(id="discard", type="text", placeholder="İstenmeyen sekans dizisi",
                              style={'marginLeft': '2px'}),
                ]
            ),

            html.Div(id="output"),
        ], style={'marginTop': "2%", 'margin': 20}
    )

    @sequence_analiz.callback(
        Output("output", "children"),
        Input("sekans", "value"),
        Input("start", "value"),
        Input("stop", "value"),
        Input("discard", "value"),
    )
    def update_output(sekans, start, stop, discard):
        sekans = str(sekans).replace(" ", "")
        sekans = str(sekans).replace("\n", "")
        sekans = str(sekans).replace("\t", "")
        sekans = str(sekans).replace("\r", "")
        sekans = str(sekans).replace("0", "")
        sekans = str(sekans).replace("1", "")
        sekans = str(sekans).replace("2", "")
        sekans = str(sekans).replace("3", "")
        sekans = str(sekans).replace("4", "")
        sekans = str(sekans).replace("5", "")
        sekans = str(sekans).replace("6", "")
        sekans = str(sekans).replace("7", "")
        sekans = str(sekans).replace("8", "")
        sekans = str(sekans).replace("9", "")
        sekans = str(sekans).upper().replace("NONE", "")

        if discard:
            sekans = sekans.replace(str(discard).upper(), "")

        if start or stop:
            sekans = sekans[start:stop]

        nuc_df = pd.DataFrame({
            'positions': [pos + 1 for pos in range(len(sekans))],
            'nucleotit': [seq for seq in sekans]
        })

        a = nuc_df[(nuc_df['nucleotit']) == "A"]
        t = nuc_df[(nuc_df['nucleotit']) == "T"]
        g = nuc_df[(nuc_df['nucleotit']) == "G"]
        c = nuc_df[(nuc_df['nucleotit']) == "C"]

        return html.Div([
            html.Hr(),

            html.P(
                f"SEKANS UZUNLUĞU: {len(sekans)}, %GC: {gc_fraction(sekans)}," +
                f"  A: {sekans.count('A')}, T: {sekans.count('T')}, G: {sekans.count('G')}, C: {sekans.count('C')}",
                style={'marginTop': '10px'}
            ),

            html.Hr(),

            html.Div([

                html.P("Grafik", className='text-primary'),

                dcc.Graph(
                    figure=ff.create_distplot(
                        [a['positions'], t['positions'], g['positions'], c['positions']],
                        ['A', 'T', 'G', 'C'],
                        bin_size=[.1, .25, .5, 1], show_hist=False,
                        curve_type="normal",
                    ).update_layout({'title': 'Nükleotit Dağılım GRAFİĞİ'})),

                dcc.Graph(figure=px.histogram(
                    nuc_df.to_dict("records"), x="nucleotit", color="nucleotit",
                    title="Nükleotit Sayıları",
                ).update_yaxes(title="Nükletotit Sayısı").update_xaxes(title="Nükleotit")),
            ])
        ], style={'marginTop': "2%", 'margin': 20})

    return HttpResponseRedirect("/laboratory/bioinformatic/app/sequence_analiz/")


def translation(request):
    external_stylesheets = ['https://cdn.jsdelivr.net/npm/bootstrap@5.3.2/dist/css/bootstrap.min.css']
    external_scripts = ['https://cdn.jsdelivr.net/npm/bootstrap@5.3.2/dist/js/bootstrap.bundle.min.js']
    translate_app = DjangoDash("translate_dna", external_stylesheets=external_stylesheets,
                               external_scripts=external_scripts, add_bootstrap_links=True, title='PROTEİN SENTEZİ')

    translate_app.layout = html.Div(

        [

            html.H4('PROTEİN SENTEZİ', className="text-primary"),

            html.A('BİYOİNFORMATİK ANASAYFA', href=HttpResponseRedirect(reverse("bioinformatic:home")).url,
                   style={'float': 'right'}),

            html.P("Sekans", style={'font-weight': 'bold'}, className="text-primary"),

            dcc.Textarea(id="sekans", placeholder="Sekans Giriniz",
                         className="form-control mb-1", style={'height': 200}),

            dcc.Input(id="start", type="number", placeholder="Sekans başlangıcı", min=0, className="mr-1"),

            dcc.Input(id="stop", type="number", placeholder="Sekans bitişi", min=0, className="mr-1"),

            dcc.Input(id="discard", type="text", placeholder="İstenmeyen sekans dizisi", className="mr-1"),

            html.P("Dönüşüm Tablosu", className="fw-bold mt-2 text-primary"),

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
            }, value="1", searchable=True, style={'marginTop': -5}),

            html.Small(" ***Dönüşüm tablosu verileri ", className='txt-bold'),

            html.A(' NCBI ', href="https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi", target="_blank",
                   style={'text-decoration': 'none'}),

            html.Small(' sitesiden alınmıştır.'),

            dcc.Checklist(id="to_stop", options={'evet': ' Stop Kodonlarını dahil et.'}),

            html.Div(id="output"),
        ], style={'margin': 30})

    @translate_app.callback(
        Output("output", "children"),
        Input("sekans", "value"),
        Input("start", "value"),
        Input("stop", "value"),
        Input("discard", "value"),
        Input("table", "value"),
        Input("to_stop", "value"),
    )
    def update_output(sekans, start, stop, discard, table, to_stop):
        sekans = str(sekans).replace(" ", "")
        sekans = str(sekans).replace("\n", "")
        sekans = str(sekans).replace("\t", "")
        sekans = str(sekans).replace("\r", "")
        sekans = str(sekans).replace("0", "")
        sekans = str(sekans).replace("1", "")
        sekans = str(sekans).replace("2", "")
        sekans = str(sekans).replace("3", "")
        sekans = str(sekans).replace("4", "")
        sekans = str(sekans).replace("5", "")
        sekans = str(sekans).replace("6", "")
        sekans = str(sekans).replace("7", "")
        sekans = str(sekans).replace("8", "")
        sekans = str(sekans).replace("9", "")
        sekans = str(sekans).upper().replace("NONE", "")

        if discard:
            sekans = sekans.replace(str(discard).upper(), "")

        if start or stop:
            sekans = sekans[start:stop]

        if to_stop:
            trans_seq = Seq(sekans).translate(table=int(table))
        else:
            trans_seq = Seq(sekans).translate(table=int(table)).replace("*", "")

        if trans_seq:

            return html.Div([

            html.Hr(),

            html.H4("SONUÇLAR"),

            html.Hr(),

            html.P(
                f" Girilen sekans uzunluğu: {len(sekans)}, %GC: {gc_fraction(sekans)}," +
                f"  A: {sekans.count('A')}, T: {sekans.count('T')}, G: {sekans.count('G')}, C: {sekans.count('C')},"
                f"  Protein Sekans Uzunluğu: {len(trans_seq)}, Stop kodonları sayısı: {trans_seq.count('*')}",
                style={'marginTop': '20px'}, className="fw-bolder"
            ),

            html.Hr(),

            html.Label("Komplement", style={'font-weight': 'bold'}),

            dcc.Textarea(
                value=str(Seq(sekans).complement()),
                className="form-control"
            ),

            html.Label("Reverse Komplement", style={'font-weight': 'bold'}),

            dcc.Textarea(
                value=str(Seq(sekans).reverse_complement()),
                className="form-control"
            ),

            html.Label("Transcribe", style={'font-weight': 'bold'}),

            dcc.Textarea(
                value=str(Seq(sekans).transcribe()),
                className="form-control"
            ),

            html.Label(f"Protein Sekansı", style={'font-weight': 'bold'}),

            dcc.Textarea(
                value=str(trans_seq),
                className="form-control"
            ),

            html.Hr(),

        ])

    return HttpResponseRedirect("/laboratory/bioinformatic/app/translate_dna/")


def Kmer_SeqSlicing(request):
    external_stylesheets = ['https://cdn.jsdelivr.net/npm/bootstrap@5.3.2/dist/css/bootstrap.min.css']
    external_scripts = ['https://cdn.jsdelivr.net/npm/bootstrap@5.3.2/dist/js/bootstrap.bundle.min.js']
    seq_input = DjangoDash('kmer_seq_slicing', external_stylesheets=external_stylesheets,
                           external_scripts=external_scripts, add_bootstrap_links=True, title='Kmer oluşturma')

    seq_input.layout = html.Div(

        [

            html.H4('Kmer oluşturma'),

            html.A('BİYOİNFORMATİK ANASAYFA', href=HttpResponseRedirect(reverse("bioinformatic:home")).url,
                   style={'float': 'right'}),

            html.P("Sekans giriniz", className='fw-bolder'),

            dcc.Textarea(id="sekans", placeholder="Sekans Giriniz",
                         style={'marginRight': '10px', 'width': '100%', 'height': '200px'}),

            dcc.Input(id="start", type="number", placeholder="Sekans başlangıcı", min=1,
                      style={'marginLeft': '2px'}),
            dcc.Input(id="stop", type="number", placeholder="Sekans bitişi", min=1,
                      style={'marginLeft': '2px'}),
            dcc.Input(id="kmer", type="number", placeholder="Kmer Uzunluğu", min=0,
                      style={'marginLeft': '2px'}),
            dcc.Input(id="discard", type="text", placeholder="İstenmeyen sekans dizisi",
                      style={'marginLeft': '2px'}),

            html.Hr(),

            html.Div(id="output"),

        ], style={'marginTop': "2%", 'margin': 20}
    )

    @seq_input.callback(
        Output("output", "children"),
        Input("sekans", "value"),
        Input("start", "value"),
        Input("stop", "value"),
        Input("kmer", "value"),
        Input("discard", "value"),
    )
    def update_output(sekans, start, stop, kmer, discard):
        sekans = str(sekans).replace(" ", "")
        sekans = str(sekans).replace("\n", "")
        sekans = str(sekans).replace("\t", "")
        sekans = str(sekans).replace("\r", "")
        sekans = str(sekans).replace("0", "")
        sekans = str(sekans).replace("1", "")
        sekans = str(sekans).replace("2", "")
        sekans = str(sekans).replace("3", "")
        sekans = str(sekans).replace("4", "")
        sekans = str(sekans).replace("5", "")
        sekans = str(sekans).replace("6", "")
        sekans = str(sekans).replace("7", "")
        sekans = str(sekans).replace("8", "")
        sekans = str(sekans).replace("9", "")
        sekans = str(sekans).upper().replace("NONE", "")

        if discard:
            sekans = sekans.replace(str(discard).upper(), "")

        if start or stop:
            sekans = sekans[start:stop]

        kmer_list = []

        def getKmers(sequence, size):
            for x in range(0, len(sequence) - size):
                yield sequence[x:x + size]

        for km in getKmers(sekans, kmer):
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

                html.P(f"SEKANS UZUNLUĞU: {len(sekans)}, %GC: {gc_fraction(sekans)}",
                       style={'marginTop': '10px'}),

                html.P("Sekanslar", className='fw-bolder'),

                dcc.Textarea(value=sekans, style={'width': '100%', 'height': '200px'}),

                html.P("Kmerler", className="fw-bolder"),

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

    return HttpResponseRedirect("/laboratory/bioinformatic/app/kmer_seq_slicing/")


def create_frame_seq(request):
    external_stylesheets = ['https://cdn.jsdelivr.net/npm/bootstrap@5.3.2/dist/css/bootstrap.min.css']
    external_scripts = ['https://cdn.jsdelivr.net/npm/bootstrap@5.3.2/dist/js/bootstrap.bundle.min.js']

    frame_seq_app = DjangoDash('create_frame_seq', external_stylesheets=external_stylesheets,
                               external_scripts=external_scripts, add_bootstrap_links=True,
                               title='Çerçeve sekans oluşturma'.upper())

    frame_seq_app.layout = html.Div(

        [

            html.H4('Çerçeve sekans oluşturma'),

            html.A('BİYOİNFORMATİK ANASAYFA', href=HttpResponseRedirect(reverse("bioinformatic:home")).url,
                   style={'float': 'right'}),

            html.P("Sekans"),

            html.Div(
                [
                    dcc.Textarea(id="sekans", placeholder="Sekans Giriniz",
                                 style={'marginRight': '10px', 'width': '100%', 'height': '200px'}),
                ]
            ),

            html.Div(
                [
                    dcc.Input(id="nuc_pos", type="number", placeholder="Nükleotit pozisyonu", min=0,
                              style={'marginLeft': '2px'}),
                    dcc.Input(id="frame_seq_len", type="number", placeholder="Çerçeve sekans uzunluğu", min=50,
                              style={'marginLeft': '2px'}),
                    dcc.Input(id="start", type="number", placeholder="Sekans başlangıcı", min=0,
                              style={'marginLeft': '2px'}),
                    dcc.Input(id="stop", type="number", placeholder="Sekans bitişi", min=0,
                              style={'marginLeft': '2px'}),

                    dcc.Input(id="discard", type="text", placeholder="İstenmeyen sekans dizisi",
                              style={'marginLeft': '2px'}),
                ]
            ),

            html.Div(id="output"),
        ], style={'marginTop': "2%"}
    )

    @frame_seq_app.callback(
        Output("output", "children"),
        Input("sekans", "value"),
        Input("nuc_pos", "value"),
        Input("frame_seq_len", "value"),
        Input("start", "value"),
        Input("stop", "value"),
        Input("discard", "value"),
    )
    def update_output(sekans, nuc_pos, frame_seq_len, start, stop, discard):
        sekans = str(sekans).replace(" ", "")
        sekans = str(sekans).replace("\n", "")
        sekans = str(sekans).replace("\t", "")
        sekans = str(sekans).replace("\r", "")
        sekans = str(sekans).replace("0", "")
        sekans = str(sekans).replace("1", "")
        sekans = str(sekans).replace("2", "")
        sekans = str(sekans).replace("3", "")
        sekans = str(sekans).replace("4", "")
        sekans = str(sekans).replace("5", "")
        sekans = str(sekans).replace("6", "")
        sekans = str(sekans).replace("7", "")
        sekans = str(sekans).replace("8", "")
        sekans = str(sekans).replace("9", "")
        sekans = str(sekans).upper().replace("NONE", "")

        if discard:
            sekans = sekans.replace(str(discard).upper(), "")

        if start or stop:
            sekans = sekans[start:stop]

        df = pd.DataFrame({
            'seq': [k for k in sekans]
        }, index=[i + 1 for i in range(len(sekans))])

        nuc_position = df.to_dict()['seq'].get(nuc_pos)

        df = pd.DataFrame({
            'seq_len': [len(sekans[:50] + str(nuc_position) + sekans[50:frame_len]) for frame_len in
                        range(50, frame_seq_len)],
            'seq': [sekans[:50] + ' ' + str(nuc_position) + ' ' + sekans[50:frame_len] for frame_len in
                    range(50, frame_seq_len)]
        })
        return html.Div([
            html.P(
                f"SEKANS UZUNLUĞU: {len(sekans)}, %GC: {gc_fraction(sekans)}, Nükleotid Pozisyonu : {nuc_position}",
                style={'marginTop': '20px'}

            ),

            html.P("Çerçeve Sekanslar", className="fw-bolder"),

            html.Div([
                dag.AgGrid(
                    id="frame_seq_table",
                    style={'width': '100%'},
                    rowData=df.to_dict("records"),
                    columnSize="sizeToFit",
                    defaultColDef={"resizable": True, "sortable": True, "filter": True, 'editable': True,
                                   "minWidth": 125},
                    dashGridOptions={'pagination': True, "rowSelection": "multiple"},
                    columnDefs=[
                        {'field': 'seq_len', 'headerName': 'UZUNLUK', 'filter': True},
                        {'field': 'seq', 'headerName': 'SEKANS', 'filter': True},
                    ]
                ),
                html.Div(id="quickstart-output"),
            ]),
            html.Hr(),
        ]),

    @frame_seq_app.callback(
        Output("quickstart-output", "children"),
        Input("frame_seq_table", "cellClicked")
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

    return HttpResponseRedirect("/laboratory/bioinformatic/app/create_frame_seq/")
