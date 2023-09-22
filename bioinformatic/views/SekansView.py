import datetime
from pathlib import Path

import Bio
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

BASE_DIR = Path(__file__).resolve().parent.parent.parent


def sequence_analiz(request):
    if request.user.is_anonymous:
        from django.conf import settings
        messages.error(request, "Lütfen Giriş Yapınız")
        return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))

    form = DNASekansForm(request.POST or None)

    if request.method == "POST":
        if form.is_valid():

            sekans = form.cleaned_data['dna'].upper()

            sekans = sekans.replace(" ", "")
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

            my_seq = Seq("{}".format(sekans))

            protein = ['M', 'R', 'N', '	D', 'E', 'Q', 'H', 'L', 'K']

            for aa in protein:
                for dna in my_seq:
                    if aa in dna:
                        return render(request, 'exception/page-404.html',
                                      {'msg': 'Sekensda protein tespit edilmiştir.'})

            seq_dash = DjangoDash('bar')

            df = pd.DataFrame(
                {
                    'nükleotit': ["ADENİN", "TİMİN", "GUANİN", "SİTOZİN"],
                    'len': [my_seq.count('A'), my_seq.count('T'), my_seq.count('G'), my_seq.count('C')],
                }
            )

            seq_dash.layout = html.Div([
                dcc.Graph(figure=px.bar(df, x="nükleotit", y="len", color="nükleotit", title="NÜKLEOTİT",
                                        labels={
                                            "nükleotit": "Nükleotit",
                                            "len": "Uzunluk"
                                        })),
            ])

            return render(request, 'bioinformatic/sekans/result.html',
                          {
                              'len': len(my_seq), 'A': my_seq.count('A'), 'G': my_seq.count('G'),
                              'C': my_seq.count('C'), 'T': my_seq.count('T'),
                              'GC': GC(my_seq).__round__(2),
                          })

        else:
            form = DNASekansForm()

    return render(request, 'bioinformatic/form.html', {'form': form, 'title': 'DNA SEKANS OKUMASI'})


def translation(request):
    if request.user.is_anonymous:
        from django.conf import settings
        messages.error(request, "Lütfen Giriş Yapınız")
        return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))

    form = TranslationForm(request.POST or None)

    if request.method == "POST":

        if form.is_valid():

            table = form.cleaned_data['table']
            sequence = form.cleaned_data['sequence'].upper()
            to_stop = form.cleaned_data['to_stop']

            sequence = sequence.replace(" ", "")
            sequence = sequence.replace("\n", "")
            sequence = sequence.replace("\t", "")
            sequence = sequence.replace("\r", "")
            sequence = sequence.replace("0", "")
            sequence = sequence.replace("1", "")
            sequence = sequence.replace("2", "")
            sequence = sequence.replace("3", "")
            sequence = sequence.replace("4", "")
            sequence = sequence.replace("5", "")
            sequence = sequence.replace("6", "")
            sequence = sequence.replace("7", "")
            sequence = sequence.replace("8", "")
            sequence = sequence.replace("9", "")

            complement = Seq(sequence).complement()
            reverse = Seq(sequence).reverse_complement()
            transcribe = Seq(sequence).transcribe()

            try:
                if to_stop is True:
                    translate = Seq(sequence).translate(table=table)
                else:
                    translate = Seq(sequence).translate(table=table).replace("*", "")

                return render(request, 'bioinformatic/translation/result.html',
                              context={
                                  'date': datetime.datetime.now(),
                                  'table': table,
                                  'len': len(sequence),
                                  'GC': GC(sequence).__round__(2),
                                  'complement': complement,
                                  'revese': reverse,
                                  'transcribe': transcribe,
                                  'translate': translate,
                              })

            except Bio.Seq.CodonTable.TranslationError:

                msg = "Codon Tablo Hatası"

                return render(request, 'exception/page-404.html', {'msg': msg})

        else:

            form = TranslationForm()

    return render(request, 'bioinformatic/form.html', {'form': form, 'title': "DNA SEKANS TRASLASYON"})


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
                    dcc.Input(id="kmer", type="number", placeholder="Kmer Uzunluğu", minLength=1,
                              style={'marginLeft': '2px'}),
                    dcc.Input(id="nuc_pos", type="number", placeholder="Nükleotit pozisyonu", minLength=1,
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
        sekans = sekans.replace(str(discard).upper(), "")
        sekans = sekans[start:stop].upper()

        def getKmers(sequence, size, step=4):
            for x in range(0, len(sequence) - size, step):
                yield sequence[x:x + size]

        df = pd.DataFrame({
            'index': [i for i in range(len(sekans))],
            'seq': [k for k in sekans]
        })

        nuc_position = df.to_dict()['seq'].get(nuc_pos)

        kmer_list = []

        if kmer:

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
