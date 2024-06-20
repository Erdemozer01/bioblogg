import base64
import os
import tempfile
from shutil import copy2
from textwrap import dedent as s
from django.shortcuts import *
import dash_daq as daq
from dash import dcc, html, Input, Output
from Bio import SeqIO
from django_plotly_dash import DjangoDash
import dash_bootstrap_components as dbc
from sklearn import datasets
from django.contrib import messages
import pandas as pd
from bioinformatic.forms import CircosForm
from bioinformatic.models import BioinformaticModel
import dash_bio as dashbio
import json
from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram

DATAPATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')

data_info = {
    os.path.join(DATAPATH, '4uft.pdb'): {
        'name': 'Measles Nucleocapsid',
        'description': dcc.Markdown(s(r'''
        The measles nucleoprotein forms a large helical complex with
        RNA... It is thought to chaperone the process of replication and
        transcription by providing a site ready for binding of the
        polymerase/phosphoprotein complex while a new RNA chain is being
        built.

        The structure includes the stable core domain of the
        nucleoprotein and a strand of RNA, but the flexible tail of
        nucleoprotein was removed in this study.
        ''')),
        'link': 'http://pdb101.rcsb.org/motm/231'
    },

    os.path.join(DATAPATH, '1yi5.pdb'): {
        'name': 'a-cobratoxin-AChBP complex',
        'description': dcc.Markdown(s(r'''

        The crystal structure of the snake long alpha-neurotoxin,
        alpha-cobratoxin, bound to the pentameric
        acetylcholine-binding protein (AChBP) from Lymnaea
        stagnalis...

        The structure unambiguously reveals the positions and
        orientations of all five three-fingered toxin molecules
        inserted at the AChBP subunit interfaces and the
        conformational changes associated with toxin binding.
        ''')),
        'link': 'https://www.rcsb.org/structure/1yi5'
    },

    os.path.join(DATAPATH, '1su4.pdb'): {
        'name': 'Calcium ATPase',
        'description': dcc.Markdown(s(r'''
        The calcium pump allows muscles to relax after... \[muscle\]
        contraction. The pump is found in the membrane of the
        sarcoplasmic reticulum. In some cases, it is so plentiful that
        it may make up 90% of the protein there. Powered by ATP, it
        pumps calcium ions back into the sarcoplasmic reticulum,
        reducing the calcium level around the actin and myosin
        filaments and allowing the muscle to relax.

        \[The structure\] has a big domain poking out on the outside
        of the sarcoplasmic reticulum, and a region that is embedded
        in the membrane, forming a tunnel to the other side.
        ''')),
        'link': 'http://pdb101.rcsb.org/motm/51'
    },

    os.path.join(DATAPATH, '1bna.pdb'): {
        'name': 'DNA',
        'description': dcc.Markdown(s(r'''
        DNA is read-only memory, archived safely inside cells. Genetic
        information is stored in an orderly manner in strands of
        DNA. DNA is composed of a long linear strand of millions of
        nucleotides, and is most often found paired with a partner
        strand. These strands wrap around each other in the familiar
        double helix...
        ''')),
        'link': 'http://pdb101.rcsb.org/motm/23'
    },

    os.path.join(DATAPATH, '1msw.pdb'): {
        'name': 'T7 RNA Polymerase',
        'description': dcc.Markdown(s(r'''
        RNA polymerase is a huge factory with many moving parts. \[The
        constituent proteins\] form a machine that surrounds DNA
        strands, unwinds them, and builds an RNA strand based on the
        information held inside the DNA. Once the enzyme gets started,
        RNA polymerase marches confidently along the DNA copying RNA
        strands thousands of nucleotides long.

        ...

        This structure includes a very small RNA polymerase that is
        made by the bacteriophage T7... \[a\] small transcription
        bubble, composed of two DNA strands and an RNA strand, is
        bound in the active site.
        ''')),
        'link': 'http://pdb101.rcsb.org/motm/40'}

}

iris_raw = datasets.load_iris()
iris = pd.DataFrame(iris_raw["data"], columns=iris_raw["feature_names"])

external_stylesheets = [dbc.themes.BOOTSTRAP]

controls = dbc.Card(
    id='control-tabs',
    className='control-tabs',
    children=[dcc.Tabs(id='tabs', value='what-is', children=[
        dcc.Tab(
            label='About',
            value='what-is',
            children=html.Div(className='control-tab', children=[
                html.H4(className='mt-3', children='What is Molecule3D?'),
                html.P('Molecule3D is a visualizer that allows you '
                       'to view biomolecules in multiple representations: '
                       'sticks, spheres, and cartoons.'),
                html.P('You can select a preloaded structure, or upload your own, '
                       'in the "Data" tab. A sample structure is also '
                       'available to download.'),
                html.P('In the "View" tab, you can change the style and '
                       'coloring of the various components of your molecule.')
            ])
        ),

        dcc.Tab(
            label='Data',
            value='upload-download',
            children=html.Div(className='control-tab', children=[
                html.Div(
                    title='download a sample data file to view',
                    children=[
                    ]
                ),

                html.Div(
                    title='Upload biomolecule to view here',

                    id='upload-field', children=[
                        dcc.Upload(
                            id='upload-data',
                            className='form-control',
                            children=html.Div([
                                'Drag and drop or click to upload a file.',
                            ]),
                            style={
                                'width': '100%',
                                'height': '70px',
                                'lineHeight': '60px',
                                'borderWidth': '1px',
                                'borderStyle': 'dashed',
                                'borderRadius': '5px',
                                'textAlign': 'center',
                                'margin': '1px',
                                'marginTop': '3%',
                                'marginBottom': '3%',
                            },
                            # Allow multiple files to be uploaded
                            multiple=True
                        ),
                        html.A(
                            html.Button(
                                "Görüntüyü İndir",
                                id="download-img",
                                className='btn btn-primary col-12 ',
                            ),
                            href=os.path.join('assets', 'sample_data', '2mru.pdb'),
                            download='2mru.pdb'
                        )
                    ]
                ),
            ])
        ),

        dcc.Tab(
            label='View',
            value='view-options',
            children=html.Div(className='control-tab', children=[
                html.Div(
                    title='Format',
                    className="app-controls-block",
                    children=[
                        html.P(
                            'Şekil Seçiniz',
                            style={
                                'font-weight': 'bold',
                                'margin-bottom': '3%',
                                'margin-top': '3%',
                            }
                        ),
                        dcc.Dropdown(
                            id='format',
                            options=[
                                {'label': 'Lineer', 'value': 'linear'},
                                {'label': 'Dairesel', 'value': 'circular'},
                            ],
                            value='linear'
                        ),
                    ],
                ),

                html.Div(
                    title='Şekil',
                    className="app-controls-block",
                    children=[
                        html.P(
                            'Şekil',
                            style={
                                'font-weight': 'bold',
                                'margin-bottom': '3%',
                                'margin-top': '3%',
                            }
                        ),
                        dcc.Dropdown(
                            id='sigil',
                            options=[
                                {'label': 'Kutu', 'value': 'BOX'},
                                {'label': 'OK', 'value': 'ARROW'},
                                {'label': 'Sekizgen', 'value': 'OCTO'},
                                {'label': 'Pürüzlü', 'value': 'JAGGY'},
                                {'label': 'Büyük OK', 'value': 'BIGARROW'},
                            ],
                            value='ARROW'
                        ),

                    ],
                ),

                html.Div(
                    title='Renk',
                    className="app-controls-block",
                    children=[
                        html.P(
                            'Renk',
                            style={
                                'font-weight': 'bold',
                                'margin-bottom': '3%',
                                'margin-top': '3%',
                            }
                        ),
                        dcc.Dropdown(
                            id='color',
                            options=[
                                {'label': 'Mavi', 'value': 'blue'},
                                {'label': 'Kırmızı', 'value': 'red'},
                                {'label': 'Yeşil', 'value': 'green'},
                                {'label': 'Mor', 'value': 'purple'},
                                {'label': 'Turuncu', 'value': 'orange'},
                                {'label': 'Koyu Yeşil', 'value': 'darkgreen'},
                            ],
                            value='blue'
                        ),

                    ],
                ),

                html.Div(
                    title='select color scheme for viewing biomolecule',
                    className="app-controls-block",
                    children=[
                        html.P(
                            'Fragment Sayısı',
                            style={
                                'font-weight': 'bold',
                                'margin-bottom': '3%',
                                'margin-top': '3%',
                            }
                        ),
                        dcc.Input(id='fragment', type='number', min=1, max=10, step=1, className='form-control',
                                  value=1),
                    ],
                ),

                html.Div(
                    className="app-controls-block",
                    children=[
                        html.P(
                            'Label',
                            style={
                                'font-weight': 'bold',
                                'margin-bottom': '3%',
                                'margin-top': '3%',
                            }
                        ),
                        daq.BooleanSwitch(id='switch', on=True),
                    ],
                ),

            ]),
        ),

    ]),
              ],
    body=True
)


def DiagramView(request):
    if request.user.is_anonymous:
        from django.conf import settings
        messages.error(request, "Lütfen Giriş Yapınız")
        return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))

    if BioinformaticModel.objects.filter(user=request.user, tool="circos").exists():
        BioinformaticModel.objects.filter(user=request.user, tool="circos").delete()

    app = DjangoDash(name='gen-diyagram-sonuc', external_stylesheets=external_stylesheets,
                     title='Genom Diyagram')

    form = CircosForm(request.POST or None, request.FILES or None)

    if request.method == "POST":
        if form.is_valid():
            in_file = form.cleaned_data["file"]

            obj = BioinformaticModel.objects.create(
                user=request.user,
                tool="diyagram",
            )

            file_obj = obj.records_files.create(file=in_file)

            try:
                records = SeqIO.read(handle=file_obj.file.path, format='genbank')

            except ValueError:

                return render(request, "exception/page-404.html", {'msg': 'Çoklu Kayıt Girdiniz'})

            app.layout = dbc.Container(
                [
                    ##NAVBAR##
                    dbc.NavbarSimple(
                        children=[
                            dbc.NavItem(dbc.NavLink("Blog", href=HttpResponseRedirect(
                                reverse("blog:anasayfa")).url, external_link=True)),
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
                        brand="Genom Diyagram",
                        brand_href=HttpResponseRedirect(reverse("bioinformatic:genome_diagram")).url,
                        color="primary",
                        dark=True,
                        brand_external_link=True
                    ),
                    html.Hr(),
                    dbc.Row(
                        [
                            dbc.Col(controls, md=4),
                            dbc.Col(dcc.Graph(id="output-graph", className="mx-auto"), md=8),
                        ],
                        align="center",
                    ),
                ],
                fluid=False, className="shadow-lg p-3 mb-5 bg-body rounded mt-3",
            )

            @app.callback(
                Output('output-graph', 'figure'),
                Input('format', 'value'),
                Input('sigil', 'value'),
                Input('color', 'value'),
                Input('switch', 'on'),
                Input('fragment', 'value'),
            )
            def GenomeDiagramUpdate(format, sigil, color, switch, fragment):

                gd_diagram = GenomeDiagram.Diagram("Yersinia pestis biovar Microtus plasmid pPCP1")
                gd_track_for_features = gd_diagram.new_track(1, name="Annotated Features")
                gd_feature_set = gd_track_for_features.new_set()

                for feature in records.features:
                    if feature.type != "gene":
                        continue
                if len(gd_feature_set) % 2 == 0:
                    color = colors.blue
                else:
                    color = colors.lightblue
                gd_feature_set.add_feature(feature, color=color, label=True)

                results = gd_diagram.draw(
                    format="linear",
                    orientation="landscape",
                    pagesize="A4",
                    fragments=4,
                    start=0,
                    end=len(records),
                )

                return

        return HttpResponseRedirect("/laboratuvarlar/bioinformatic-laboratuvari/app/gen-diyagram-sonuc/")

    return render(request, "bioinformatic/form.html", {'form': form, 'title': 'Genom Diyagram'})
