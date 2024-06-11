import mimetypes
import os
from django.http.response import HttpResponse, HttpResponseRedirect
from pathlib import Path
from bioinformatic.models.bioinformatic import BioinformaticModel
from dash import html, Output, Input, dcc
from django_plotly_dash import DjangoDash
import dash_bootstrap_components as dbc


BASE_DIR = Path(__file__).resolve().parent.parent.parent


def EntrezDownload(request):
    external_stylesheets = ['https://cdn.jsdelivr.net/npm/bootstrap@5.3.2/dist/css/bootstrap.min.css']
    external_scripts = ['https://cdn.jsdelivr.net/npm/bootstrap@5.3.2/dist/js/bootstrap.bundle.min.js']

    app = DjangoDash('entrez', external_stylesheets=external_stylesheets, external_scripts=external_scripts,
                     title='Genbank veya Fasta Dosyası indir')

    controls = dbc.Card(
        [

            html.Div(
                [
                    dbc.Label("Email"),
                    dbc.Input(placeholder="Email", type="email", id="email"),
                ]
            ),
            html.Div(
                [
                    dbc.Label("VERİ TABANI"),
                    dcc.Dropdown(
                        id="select",
                        options={
                            "Protein": "Protein",
                            "Nucleotide": "Nucleotide",
                        },
                    ),
                ]
            ),
            html.Div(
                [
                    dbc.Label("İD NUMARASI"),
                    dbc.Input(placeholder="EU490707, 6273291", type="text", id="id"),
                ]
            ),
            html.Div(
                [
                    dbc.Label("Tris"),
                    dbc.Input(placeholder="Tris", type="number", id="Tris", value=0, min=0),
                ]
            ),
            html.Div(
                [
                    dbc.Label("Magnezyum"),
                    dbc.Input(placeholder="Magnezyum", type="number", id="Mg", value=0, min=0),
                ]
            ),
            html.Div(
                [
                    dbc.Label("dNTPs"),
                    dbc.Input(placeholder="dNTPs", type="number", id="dNTPs", value=0, min=0),
                ]
            ),
            html.Div(
                [
                    dbc.Label("Tuz"),
                    dbc.Input(placeholder="Tuz", type="number", id="saltcorr", value=5, min=0, max=6),
                ]
            ),
        ],
        body=True,
    )

    app.layout = dbc.Container(
        [
            html.H3("Primer Erime sıcaklığı"),
            html.Hr(),
            dbc.Row(
                [
                    dbc.Col(controls, md=6),
                    dbc.Col(html.Div(id="output"), md=4),
                ],
                align="center",
            ),
        ],
        className="m-4"
    )

    @app.callback(
        Output("output", 'children'),
        Input("sekans", "value"),
        Input("Na", "value"),
        Input("K", "value"),
        Input("Tris", "value"),
        Input("Mg", "value"),
        Input("dNTPs", "value"),
        Input("saltcorr", "value")
    )
    def temp_melting(sekans, Na, K, Tris, Mg, dNTPs, saltcorr):
        if sekans is not None:
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

            primer_len = len(sekans)

            temperature_melting = mt.Tm_NN(seq=sekans, Na=Na, K=K, Tris=Tris, Mg=Mg, dNTPs=dNTPs, saltcorr=saltcorr)

            return html.Div([
                html.P(f"Sekans Uzunluğu : {primer_len}"),
                html.P(f"Erime Sıcaklığı : {temperature_melting}"),
            ])

        else:
            return "Sekans girilmedi"

    return HttpResponseRedirect("laboratuvarlar/bioinformatic-laboratuvari/app/temp-melt/")


def download_file(request):
    obj = BioinformaticModel.objects.filter(user=request.user)
    format = [j.writing_file_format for j in obj][0]
    # Define text file name
    filename = f'{request.user}_{format}_file.{format}'
    # Define the full file path
    filepath = os.path.join(BASE_DIR, 'media', f'{request.user}_{format}_file.{format}')
    # Open the file for reading content
    path = open(filepath, 'r')
    # Set the mime type
    mime_type, _ = mimetypes.guess_type(filepath)
    # Set the return value of the HttpResponse
    response = HttpResponse(path, content_type=mime_type)
    # Set the HTTP header for sending to browser
    response['Content-Disposition'] = "attachment; filename=%s" % filename
    try:
        # Return the response value
        return response
    finally:
        path.close()
        os.remove(filepath)
        obj.delete()
