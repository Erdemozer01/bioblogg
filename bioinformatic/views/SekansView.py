import datetime
from pathlib import Path

import Bio
import dash_bootstrap_components.themes
import pandas as pd
import plotly.express as px
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from dash import html, dcc, Input, Output
from django.contrib import messages
from django.shortcuts import *
from django_plotly_dash import DjangoDash
from bioinformatic.forms import DNASekansForm, TranslationForm
import dash_ag_grid as dag
from dash_bootstrap_components._components import Checkbox

BASE_DIR = Path(__file__).resolve().parent.parent.parent


def sequence_analiz(request):
    external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
    external_scripts = ['https://cdn.jsdelivr.net/npm/bootstrap@5.3.2/dist/js/bootstrap.bundle.min.js']
    sequence_analiz = DjangoDash("sequence_analiz", external_stylesheets=external_stylesheets,
                                 external_scripts=external_scripts)

    sequence_analiz.layout = html.Div(

        [

            html.H4('Sekans Analiz'),

            html.A('BİYOİNFORMATİK ANASAYFA', href=HttpResponseRedirect(reverse("bioinformatic:home")).url,
                   style={'float': 'right'}),

            html.P("Sekans"),

            html.Div(
                [
                    dcc.Textarea(id="sekans", placeholder="Sekans Giriniz",
                                 style={'marginRight': '10px', 'width': '100%', 'height': '200px'}),
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

    @sequence_analiz.callback(
        Output("output", "children"),
        Input("sekans", "value"),
        Input("start", "value"),
        Input("stop", "value"),
        Input("discard", "value"),
    )
    def update_output(sekans, start, stop, discard):
        sekans = sekans.replace(' ', '')
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
        sekans = sekans.upper()

        if discard:
            sekans = sekans.replace(str(discard).upper(), "")

        if start or stop:
            sekans = sekans[start:stop].upper()

        nuc_df = pd.DataFrame({
            'positions': [pos for pos in range(len(sekans))],
            'seq': [seq for seq in sekans]
        })

        return html.Div([
            html.Hr(),
            html.P(
                f"SEKANS UZUNLUĞU: {len(sekans)}, %GC: {GC(sekans).__round__(2)}," +
                f"  A: {sekans.count('A')}, T: {sekans.count('T')}, G: {sekans.count('G')}, C: {sekans.count('C')}",
                style={'marginTop': '20px'}

            ),

            html.Hr(),

            html.Label("Sekanslar"),

            html.Div([
                dcc.Textarea(
                    value=sekans,
                    style={'width': '100%', 'height': '200px'}
                )
            ]),

            html.Hr(),

            html.Div([
                html.Label("Grafik"),
                dcc.Graph(figure=px.histogram(nuc_df, x="seq", y="positions", color="seq")),
            ])
        ])

    return HttpResponseRedirect("/laboratory/bioinformatic/app/sequence_analiz/")


def translation(request):
    external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
    external_scripts = ['https://cdn.jsdelivr.net/npm/bootstrap@5.3.2/dist/js/bootstrap.bundle.min.js']
    translate_app = DjangoDash("translate_dna", external_stylesheets=external_stylesheets,
                               external_scripts=external_scripts)

    translate_app.layout = html.Div(

        [

            html.H4('PROTEİN SENTEZİ'),

            html.A('BİYOİNFORMATİK ANASAYFA', href=HttpResponseRedirect(reverse("bioinformatic:home")).url,
                   style={'float': 'right'}),

            html.Div(
                [

                    html.P("Sekans", style={'font-weight': 'bold'}),

                    dcc.Textarea(id="sekans", placeholder="Sekans Giriniz",
                                 style={'marginRight': '10px', 'width': '100%', 'height': '200px'}),

                    dcc.Input(id="start", type="number", placeholder="Sekans başlangıcı", min=0,
                              style={'marginLeft': '2px'}),

                    dcc.Input(id="stop", type="number", placeholder="Sekans bitişi", min=0,
                              style={'marginLeft': '2px'}),

                    dcc.Input(id="discard", type="text", placeholder="İstenmeyen sekans dizisi",
                              style={'marginLeft': '2px'}),

                    html.Hr(),

                    html.P("Dönüşüm Tablosu", style={'font-weight': 'bold'}),

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
                    html.A(' NCBI ', href="https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi", target="_blank",
                           style={'text-decoration': 'none'}),
                    html.Small(' sitesiden alınmıştır.'),

                    dcc.Checklist(id="to_stop", options={'evet': 'Stop Kodonlarını dahil et.'}),

                    html.Br(),

                ]
            ),

            html.Div(id="output"),
        ], style={'marginTop': "2%"}
    )

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
        sekans = sekans.replace(' ', '')
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
        sekans = sekans.upper()

        if discard:
            sekans = sekans.replace(str(discard).upper(), "")

        if start or stop:
            sekans = sekans[start:stop].upper()

        if to_stop:
            trans_seq = Seq(sekans).translate(table=int(table))
        else:
            trans_seq = Seq(sekans).translate(table=int(table)).replace("*", "")

        return html.Div([
            html.Hr(),
            html.P(
                f" Girilen sekans uzunluğu: {len(sekans)}, %GC: {GC(sekans).__round__(2)}," +
                f"  A: {sekans.count('A')}, T: {sekans.count('T')}, G: {sekans.count('G')}, C: {sekans.count('C')},"
                f"  Protein Sekans Uzunluğu: {len(trans_seq)}, Stop kodonları sayısı: {trans_seq.count('*')}",

                style={'marginTop': '20px'}

            ),

            html.Hr(),

            html.Label("Protein Sekansı", style={'font-weight': 'bold'}),

            html.Div([
                dcc.Textarea(
                    value=str(trans_seq),
                    style={'width': '100%', 'height': '200px'}
                )
            ]),

            html.Hr(),
        ])

    return HttpResponseRedirect("/laboratory/bioinformatic/app/translate_dna/")


def Kmer_SeqSlicing(request):
    external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
    seq_input = DjangoDash('kmer_seq_slicing', external_stylesheets=external_stylesheets)

    seq_input.layout = html.Div(

        [

            html.H4('Kmer ve Sekans Kesme'),

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
                    dcc.Input(id="start", type="number", placeholder="Sekans başlangıcı", min=0,
                              style={'marginLeft': '2px'}),
                    dcc.Input(id="stop", type="number", placeholder="Sekans bitişi", min=0,
                              style={'marginLeft': '2px'}),
                    dcc.Input(id="kmer", type="number", placeholder="Kmer Uzunluğu", min=0,
                              style={'marginLeft': '2px'}),
                    dcc.Input(id="nuc_pos", type="number", placeholder="Nükleotit pozisyonu", min=0,
                              style={'marginLeft': '2px'}),
                    dcc.Input(id="discard", type="text", placeholder="İstenmeyen sekans dizisi",
                              style={'marginLeft': '2px'}),
                ]
            ),

            html.Div(id="output"),
        ], style={'marginTop': "2%"}
    )

    @seq_input.callback(
        Output("output", "children"),
        Input("sekans", "value"),
        Input("start", "value"),
        Input("stop", "value"),
        Input("kmer", "value"),
        Input("nuc_pos", "value"),
        Input("discard", "value"),
    )
    def update_output(sekans, start, stop, kmer, nuc_pos, discard):
        sekans = sekans.replace(' ', '')
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

        if discard:
            sekans = sekans.replace(str(discard).upper(), "")

        if start or stop:
            sekans = sekans[start:stop].upper()

        df = pd.DataFrame({
            'index': [i for i in range(len(sekans))],
            'seq': [k for k in sekans]
        })

        nuc_position = df.to_dict()['seq'].get(nuc_pos)

        kmer_list = []

        if kmer:
            def getKmers(sequence, size, step=4):
                for x in range(0, len(sequence) - size, step):
                    yield sequence[x:x + size]

            for km in getKmers(sekans, kmer):
                if km.count(str(nuc_position)) > 0:
                    kmer_list.append(km)

        df_kmer = pd.DataFrame(
            {
                f'{nuc_position} Sayısı': [i.count(nuc_position) for i in kmer_list],
                'kmer': kmer_list,
                '%gc': [GC(gc).__round__(2) for gc in kmer_list]
            }
        )

        return html.Div([
            html.P(
                f"SEKANS UZUNLUĞU: {len(sekans)}, %GC: {GC(sekans).__round__(2)}, Nükleotid Pozisyonu : {nuc_position}",
                style={'marginTop': '20px'}

            ),

            html.Label("Sekanslar"),

            html.Div([
                dcc.Textarea(
                    value=sekans,
                    style={'width': '100%', 'height': '200px'}
                )
            ]),

            html.Label("Kmerler", style={'text': 'bold'}),

            html.Div(
                [
                    dag.AgGrid(
                        style={'width': '100%'},
                        rowData=df_kmer.to_dict("records"),
                        dashGridOptions={'pagination': True},
                        columnDefs=[
                            {'field': f'{nuc_position} Sayısı'},
                            {'field': 'kmer', 'headerName': 'KMERS', 'filter': True},
                            {'field': '%gc', 'headerName': '%GC', 'filter': True},
                        ]
                    ),
                ],
            ),
            html.Label("Grafik"),
            html.Div([
                dcc.Graph(figure=px.scatter(df_kmer, x=df_kmer[f'{nuc_position} Sayısı'], color=df_kmer['kmer']))
            ])
        ])

    return HttpResponseRedirect("/laboratory/bioinformatic/app/kmer_seq_slicing/")
