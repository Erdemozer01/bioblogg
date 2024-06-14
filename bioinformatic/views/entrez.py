##### MEHMET ERDEM ÖZER, mozer232@posta.pau.edu.tr ######
from Bio import Entrez
from bioinformatic.forms import ArticleForm, EntrezSelectForm
from django.shortcuts import *
from django_plotly_dash import DjangoDash
from dash import dcc, html, Input, Output
import dash_bootstrap_components as dbc
import pandas as pd
from Bio import Medline


def entrez_tools(request):
    form = EntrezSelectForm(request.POST or None)

    if request.method == "POST":

        if form.is_valid():
            select = form.cleaned_data["select"]
            entrez_email = form.cleaned_data['email']
            Entrez.email = entrez_email

            if select == 'art':

                external_stylesheets = [dbc.themes.MORPH]

                app = DjangoDash('entrez-makale', external_stylesheets=external_stylesheets,
                                 title='MAKALE')

                app.layout = html.Div([

                    ## NAVBAR ##

                    dbc.NavbarSimple(
                        children=[
                            dbc.NavItem(dbc.NavLink("Blog", href=HttpResponseRedirect(
                                reverse("bioinformatic:entrez_tools")).url, external_link=True)),
                            dbc.DropdownMenu(
                                children=[
                                    dbc.DropdownMenuItem("Biyoinformatik",
                                                         href=HttpResponseRedirect(reverse("bioinformatic:home")).url,
                                                         external_link=True),
                                    dbc.DropdownMenuItem("Biyoinformatik",
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
                        brand="Entrez Araçları",
                        brand_href=HttpResponseRedirect(reverse("bioinformatic:entrez_tools")).url,
                        color="primary",
                        dark=True,
                        brand_external_link=True
                    ),

                    html.Hr(),

                    html.Div(
                        [
                            html.P('GÜNCEL MAKALE ARAMA', className='text-primary', style={'marginLeft': "1%"}),

                            html.Hr(),

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

                )
                def display_output(term):

                    pub_date = []

                    if term:

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
                                "BAŞLIK": [html.A(record.get("TI"),
                                                  href=f'https://pubmed.ncbi.nlm.nih.gov/{record.get("PMID")}',
                                                  target="_blank")
                                           for record in records],
                                "Yayın Tarihi": pub_date
                            }
                        )

                        return [dbc.Table.from_dataframe(df, striped=True, bordered=True, hover=True)]

                    else:
                        return html.P("ARANACAK TERİM GİRİNİZ", className="text-center text-danger")

                return HttpResponseRedirect("/laboratuvarlar/bioinformatic-laboratuvari/app/entrez-makale/")

            elif select == 'nuc':

                external_stylesheets = [dbc.themes.MORPH]

                app = DjangoDash('entrez-makale', external_stylesheets=external_stylesheets,
                                 title='MAKALE')

                app.layout = html.Div([

                    dbc.NavbarSimple(
                        children=[
                            dbc.NavItem(dbc.NavLink("Entrez Araçları", href=HttpResponseRedirect(
                                reverse("bioinformatic:entrez_tools")).url, external_link=True)),
                            dbc.DropdownMenu(
                                children=[

                                    dbc.DropdownMenuItem("BLOG",
                                                         href=HttpResponseRedirect(reverse("blog:anasayfa")).url,
                                                         external_link=True),
                                    dbc.DropdownMenuItem("Laboratuvarlar",
                                                         href=HttpResponseRedirect(reverse("lab_home")).url,
                                                         external_link=True),
                                ],
                                nav=True,
                                in_navbar=True,
                                label="DİĞER",
                            ),
                        ],
                        brand="Biyoinformatik",
                        brand_href=HttpResponseRedirect(reverse("bioinformatic:home")).url,
                        color="primary",
                        dark=True,
                        brand_external_link=True
                    ),

                    html.Hr(),

                    html.Div(
                        [
                            html.H4('20 GÜNCEL MAKALE ARAMA', className='text-primary'),

                            html.Hr(),

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

    return render(request, "bioinformatic/form.html", {'form': form, 'title': "ENTREZ ARAÇLARI"})
