##### MEHMET ERDEM ÖZER, mozer232@posta.pau.edu.tr ######
from Bio import Entrez
from django.shortcuts import *
from django_plotly_dash import DjangoDash
from dash import dcc, html, Input, Output
import dash_bootstrap_components as dbc
import pandas as pd
from Bio import Medline


def ArticleView(request):

    external_stylesheets = [dbc.themes.BOOTSTRAP]

    app = DjangoDash('makale-arama', external_stylesheets=external_stylesheets,
                     title='PubMed GÜNCEL MAKALE ARAMA', update_title="Makale Aranıyor...", add_bootstrap_links=True)

    app.layout = html.Div([

        ## NAVBAR ##

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
                                             href=HttpResponseRedirect(reverse("biyoistatislik")).url,
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
            brand="PubMed GÜNCEL MAKALE ARAMA",
            brand_href=HttpResponseRedirect(reverse("bioinformatic:article_search")).url,
            color="primary",
            dark=True,
            brand_external_link=True, sticky='top',
            className="shadow-lg p-3 bg-body rounded"
        ),

        html.Div(
            [
                html.Label("Email", className='mb-2 mt-5'),

                dcc.Input(id="email", type="email", placeholder="Email adresi",
                          className="form-control"),
                html.P("*** PubMed veritabanı için gerekli *** ", className="text-danger text-small"),

                html.Label("ARADIĞINIZ TERİM", className='mb-2 mt-3'),

                dcc.Input(id="term", type="text", placeholder="ARADIĞINIZ TERİM",
                          className="form-control"),

            ]
        ),

        html.Hr(),

        html.Div(id="output"),
    ],

        style={'marginTop': "2%"}, className='float-center m-5'
    )

    @app.callback(
        Output('output', 'children'),
        Input('term', 'value'),
        Input('email', 'value'),
    )
    def display_output(term, email):

        pub_date = []

        if term and email:

            Entrez.email = email

            stream = Entrez.esearch(db="pubmed", term=term)

            id_lists = Entrez.read(stream)

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
                    "Yayın Tarihi": pub_date
                }
            )

            return [dbc.Table.from_dataframe(df, striped=True, bordered=True, hover=True)]

        else:
            return html.P("ARANACAK TERİM GİRİNİZ", className="text-center text-danger")

    return HttpResponseRedirect("/laboratuvarlar/bioinformatic-laboratuvari/app/makale-arama/")
