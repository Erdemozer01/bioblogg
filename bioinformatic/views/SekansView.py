import datetime
from pathlib import Path

import Bio
import pandas as pd
import plotly.express as px
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from dash import html, dcc
from django.contrib import messages
from django.shortcuts import render, redirect
from django_plotly_dash import DjangoDash
from bioinformatic.forms import DNASekansForm, TranslationForm, SequenceSlicingForm

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


def SequenceSlicing(request):

    form = SequenceSlicingForm(request.POST or None)

    if request.method == "POST":

        if form.is_valid():

            start = form.cleaned_data['start']
            stop = form.cleaned_data['finish']

            sekans = form.cleaned_data['seq'].upper()

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

            my_seq = Seq(sekans)

            result = my_seq[start:stop + 1]

            df = pd.DataFrame({
                'Sekans': [seq for seq in sekans]
            })

            pd.set_option('display.max_rows', None)

            app = DjangoDash("seq_slice")

            app.layout = html.Div(
                [
                    html.P('Sekans Pozisyonları'),
                    dcc.Textarea(
                        id='textarea-example',
                        value=f'{df}',
                        style={'width': '100%', 'height': 300},
                    ),
                ]
            )

            return render(request, "bioinformatic/sekans/slice.html", {'result': result, 'start': start, 'stop': stop})

        else:

            form = SequenceSlicingForm()

    return render(request, "bioinformatic/form.html", {'form': form})
