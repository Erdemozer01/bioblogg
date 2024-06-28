##### MEHMET ERDEM ÖZER, mozer232@posta.pau.edu.tr ######
from Bio import Entrez
from django.shortcuts import *
from django_plotly_dash import DjangoDash
from dash import dcc, html, Input, Output
import dash_bootstrap_components as dbc
import pandas as pd
from Bio import Medline, SeqIO


def EntrezToolsView(request):
    external_stylesheets = [dbc.themes.BOOTSTRAP]

    app = DjangoDash('entrez-tools', external_stylesheets=external_stylesheets,
                     title='Entrez araçları', update_title="Aranıyor...", add_bootstrap_links=True)

    app.layout = dbc.Card(
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
                                                     reverse("biyoistatislik")).url,
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
                brand_href=HttpResponseRedirect(reverse("bioinformatic:entrez_tools")).url,
                color="primary",
                dark=True,
                brand_external_link=True,
                sticky='top',
                className="shadow-lg bg-body rounded mt-1 mb-1 mr-1 ml-1",
            ),

            dbc.Card(
                [

                    html.Label("Email adresi giriniz", className='fw-bolder mb-1 mt-1'),
                    dcc.Input(id="email", type="email", placeholder="Email adresi",
                              className="form-control"),

                    html.Label("Arama türü seçiniz", className="fw-bolder mt-2"),
                    dcc.Dropdown(
                        id='search-type',
                        options=[
                            {'label': 'MAKALE', 'value': 'article'},
                            {'label': 'Nükleotit (GENBANK)', 'value': 'gb_nuc'},
                        ],
                        value='article',
                    ),

                    html.P("*** PubMed veritabanı için gerekli ***",
                           className="text-danger text-small"),

                    html.Label("Aramak istediğiniz kelime yada terim giriniz", className='fw-bolder mb-2 mt-1'),

                    dcc.Input(id="term", type="text", placeholder="ARADIĞINIZ TERİM",
                              className="form-control"),

                    html.Label("Bulunacak sonuç sayısı", className='fw-bolder mb-2 mt-1'),

                    dcc.Input(id="retmax", type="number", placeholder="Bulunacak sonuç sayısı",
                              className="form-control", value=20, min=20),

                    html.Hr(),

                ], className="shadow-lg p-3 bg-body rounded mr-1 ml-1 mb-2"
            ),

            dbc.Card(id="output", children=[], className="shadow-lg p-3 bg-body rounded mr-1 ml-1 mb-2"),
        ],
    )

    @app.callback(
        Output('output', 'children'),
        Input('search-type', 'value'),
        Input('term', 'value'),
        Input('email', 'value'),
        Input('retmax', 'value'),
    )
    def display_output(type, term, email, retmax):

        if term is None and email is not None:
            return html.P("Lütfen aramak istediğiniz terimi giriniz.", className="text-center text-danger")

        elif term and email is not None:

            Entrez.email = email

            if type == 'article':

                pub_date = []

                search = Entrez.esearch(db="pubmed", term=term, retmax=retmax)

                id_lists = Entrez.read(search)

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
                        "İndex": [i + 1 for i in range(len(pub_date))],
                        "Makale Başlıkları": [html.A(record.get("TI"),
                                                     href=f'https://pubmed.ncbi.nlm.nih.gov/{record.get("PMID")}',
                                                     target="_blank", style={'text-decoration': 'none'})
                                              for record in records],
                        "Tarih": pub_date
                    }
                )

                return dbc.Table.from_dataframe(df, striped=True, bordered=True, hover=True)

            elif type == 'gb_nuc':

                stream = Entrez.esearch(db="nuccore", term=term, retmax=retmax)
                record = Entrez.read(stream)
                gi_list = record["IdList"]
                gi_str = ",".join(gi_list)
                stream = Entrez.efetch(db="nuccore", id=gi_str, rettype="gb", retmode="text")
                records = SeqIO.parse(stream, "gb")

                df = pd.DataFrame(
                    [
                        {
                            'İndex': [i + 1 for i in range(len([i.name]))],
                            'Genbank İD': i.name,
                            'BAŞLIK': html.A(
                                children=[i.description],
                                href=f'https://www.ncbi.nlm.nih.gov/nuccore/{i.name}',
                                target="_blank", style={'text-decoration': 'none'}),
                            'TARİH': i.annotations["date"]
                        }

                        for i in records
                    ]
                )

                return dbc.Table.from_dataframe(df, striped=True, bordered=True, hover=True)

        else:

            return html.P("Email adresi ve terim girmediniz", className="text-center text-danger")

    return HttpResponseRedirect("/laboratuvarlar/bioinformatic-laboratuvari/app/entrez-tools/")
