import datetime
import dash_bio as dashbio

import pandas as pd
import plotly.express as px
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from dash import html, dash_table, dcc
from django.shortcuts import render
from django_plotly_dash import DjangoDash
from bioinformatic.forms.dna import DNASekansForm


def sequence_analiz(request):
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

            protein = ['M', 'R', 'N', '	D', 'E', 'Q', 'H']

            for aa in protein:
                for dna in my_seq:
                    if aa in dna:
                        return render(request, 'exception/page-404.html',
                                      {'msg': 'Sekensda protein tespit edilmiştir.'})

            app = DjangoDash('sekans')

            df = pd.DataFrame(
                {
                    'NÜKLEOTİT': ["ADENİN", "TİMİN", "GUANİN", "SİTOZİN"],
                    'UZUNLUK': [my_seq.count('A'), my_seq.count('T'), my_seq.count('G'), my_seq.count('C')]
                }
            )

            app.layout = html.Div([
                html.Div(children='Grafik'),
                dcc.Graph(figure=px.bar(df, x="NÜKLEOTİT", y="UZUNLUK", color="NÜKLEOTİT")),
            ])

            return render(request, 'bioinformatic/sekans/result.html',
                          {
                              'len': len(my_seq), 'A': my_seq.count('A'), 'G': my_seq.count('G'),
                              'C': my_seq.count('C'), 'T': my_seq.count('T'),
                              'GC': GC(my_seq).__round__(2),
                          })

        else:
            form = DNASekansForm()

    return render(request, 'bioinformatic/sekans/analiz.html', {'form': form})
