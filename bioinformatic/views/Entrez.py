##### MEHMET ERDEM ÖZER, mozer232@posta.pau.edu.tr ######
from Bio import Entrez
from django.shortcuts import *
from django_plotly_dash import DjangoDash
from dash import dcc, html, Input, Output
import dash_bootstrap_components as dbc
import pandas as pd
from Bio import Medline, SeqIO
from django.views.generic import *


class EntrezView(TemplateView):
    template_name = 'bioinformatic/entrez.html'

    def get(self, request, *args, **kwargs):

        external_stylesheets = [dbc.themes.BOOTSTRAP]

        app = DjangoDash(
            name='entrez',
            external_stylesheets=external_stylesheets,
            title='Entrez araçları',
            update_title="Aranıyor...",
            add_bootstrap_links=True
        )

        app.layout = dbc.Card(
            [
                dbc.CardBody(
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

                                ),
                            ],
                            brand="Entrez Araçları",
                            brand_href=HttpResponseRedirect(reverse("bioinformatic:entrez_view")).url,
                            color="primary",
                            dark=True,
                            brand_external_link=True,
                            sticky='top',
                            className="shadow-lg bg-body rounded mt-1 mb-1 mr-1 ml-1",
                        ),

                        dbc.Card(
                            [

                                html.Label("Arama türü seçiniz", className="fw-bolder mt-1"),
                                dcc.Dropdown(
                                    id='search-type',
                                    options=[
                                        {'label': 'MAKALE', 'value': 'article'},
                                        {'label': 'Nükleotit (GENBANK)', 'value': 'gb_nuc'},
                                    ],
                                    value='article',
                                ),

                                html.Label("Aramak istediğiniz kelime yada terim giriniz", className='fw-bolder mb-2'),

                                dcc.Input(id="term", type="text", placeholder="ARADIĞINIZ TERİM",
                                          className="form-control"),

                                html.Label("İstenilen Sonuç Sayısı", className='fw-bolder mb-2 mt-1'),

                                dcc.Input(id="retmax", type="number", placeholder="İstenilen Sonuç Sayısı",
                                          className="form-control", value=20, min=20),

                                html.Div(id="count", children=[], className="fw-bolder mt-1"),

                                html.Hr(),

                            ], className="shadow-lg p-3 bg-body rounded mr-1 ml-1 mb-2"
                        ),

                        dbc.Card(id="output", children=[], className="shadow-lg p-3 bg-body rounded mr-1 ml-1 mb-2"),
                    ],
                )
            ]
        )

        @app.callback(
            Output('output', 'children'),
            Output('count', 'children'),

            Input('search-type', 'value'),
            Input('term', 'value'),
            Input('retmax', 'value')
        )
        def display_output(type, term, retmax):

            if term is None:
                return html.P("Lütfen aramak istediğiniz terimi giriniz.",
                              className="fw-bolder text-center text-danger"), ""

            elif term:

                Entrez.email = "A.N.Other@example.com"

                if type == 'article':

                    pub_date = []

                    search = Entrez.esearch(db="pubmed", term=term, retmax=retmax)

                    id_lists = Entrez.read(search)

                    count = f'Toplam {id_lists["Count"]} sonuçtan {retmax} tanesi bulunmuştur.'

                    idlist = id_lists["IdList"]

                    stream_fetch = Entrez.efetch(db="pubmed", id=",".join(idlist), rettype="medline",
                                                 retmode="text")

                    records = Medline.parse(stream_fetch)

                    records = list(records)

                    for i in records:
                        if i.get("DP"):
                            pub_date.append(i.get("DP"))
                        else:
                            pub_date.append("-")

                    df = pd.DataFrame(
                        {
                            "Makale Başlıkları": [html.A(record.get("TI"),
                                                         href=f'https://pubmed.ncbi.nlm.nih.gov/{record.get("PMID")}',
                                                         target="_blank", style={'text-decoration': 'none'})
                                                  for record in records],
                            "Tarih": pub_date
                        }
                    )

                    return dbc.Table.from_dataframe(df, bordered=True, hover=True, index=True,
                                                    responsive=True), count

                elif type == 'gb_nuc':

                    stream = Entrez.esearch(db="nuccore", term=term, retmax=retmax)
                    record = Entrez.read(stream)
                    gi_list = record["IdList"]
                    gi_str = ",".join(gi_list)
                    stream = Entrez.efetch(db="nuccore", id=gi_str, rettype="gb", retmode="text")
                    records = SeqIO.parse(stream, "gb")
                    count = f'Toplam {record["Count"]} sonuçtan {retmax} tanesi bulunmuştur.'

                    df_nuc = pd.DataFrame(
                        [
                            {
                                'Genbank İD': rec.name,

                                'Başlık': html.A(
                                    children=[rec.description],
                                    href=f'https://www.ncbi.nlm.nih.gov/nuccore/{rec.name}',
                                    target="_blank", style={'text-decoration': 'none'}),

                                'Tarih': rec.annotations["date"]
                            }

                            for rec in records
                        ]
                    )

                    return dbc.Table.from_dataframe(df_nuc, bordered=True, hover=True, index=True,
                                                    responsive=True), count

        return super().get(request, *args, **kwargs)

