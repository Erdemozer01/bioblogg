import parmed, os
from django.shortcuts import *
from django.contrib import messages
from bioinformatic.forms import SingleMoleculeViewForm, MultiMoleculeViewForm
from bioinformatic.models import BioinformaticModel
from django_plotly_dash import DjangoDash
from dash import dcc, html, dash_table, Input, Output
import dash_bootstrap_components as dbc
import dash_bio
from dash_bio.utils import PdbParser, create_mol3d_style, ngl_parser
import pandas as pd
from Bio.PDB import parse_pdb_header
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent.parent


def single_molecule_view(request):
    if request.user.is_anonymous:
        from django.conf import settings
        messages.error(request, "Lütfen Giriş Yapınız")
        return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))

    if BioinformaticModel.objects.filter(user=request.user, tool="molecule_view").exists():
        BioinformaticModel.objects.filter(user=request.user, tool="molecule_view").delete()

    external_stylesheets = [dbc.themes.BOOTSTRAP]

    app = DjangoDash('molecule-3d-viewer', external_stylesheets=external_stylesheets,
                     title='3D MOLEKÜL GÖRÜNTÜLEME', add_bootstrap_links=True)

    form = SingleMoleculeViewForm(request.POST or None, request.FILES or None)

    if request.method == "POST":

        if form.is_valid():

            if request.FILES['file'].size > 9 * 1024 * 1024:
                messages.error(request, "Dosya boyutu en fazla 9mb olmalıdır.")
                return HttpResponseRedirect(request.path)

            file = form.cleaned_data["file"]

            obj = BioinformaticModel.objects.create(
                user=request.user,
                tool="molecule_view",
            )

            file_obj = obj.records_files.create(file=file)

            handle = open(file_obj.file.path, 'r', encoding='utf-8')

            structure = parse_pdb_header(handle)

            try:

                parser = PdbParser(file_obj.file.path)

                data = parser.mol3d_data()

                styles = create_mol3d_style(
                    data['atoms'], visualization_type='cartoon', color_element='residue'
                )

                df = pd.DataFrame(data["atoms"])

                df['positions'] = df['positions'].apply(lambda x: ', '.join(map(str, x)))

            except ValueError:
                messages.error(request, "Sayfayı yenilediğiniz İçin veriler kaybolmuştur")

                return redirect('bioinformatic:single_molecule_3d_view')

            except AttributeError:
                messages.error(request, "Beklenmedik hata oluştu")

                return redirect('bioinformatic:single_molecule_3d_view')

            except parmed.exceptions.MoleculeError:
                messages.error(request, "Beklenmedik hata oluştu")

                return redirect('bioinformatic:single_molecule_3d_view')

            except parmed.exceptions.FormatNotFound:

                messages.error(request, "Hatalı Dosya Formatı")

                return redirect('bioinformatic:single_molecule_3d_view')

            columns = [
                {'name': 'Seri', 'id': 'serial'},
                {'name': 'Adı', 'id': 'name'},
                {'name': 'ELEMENT', 'id': 'elem'},
                {'name': 'Pozisyon', 'id': 'positions'},
                {'name': 'Kütle Büyüklüğü', 'id': 'mass_magnitude'},
                {'name': 'İndex', 'id': 'residue_index'},
                {'name': 'Bölge Adı', 'id': 'residue_name'},
                {'name': 'Zincir', 'id': 'chain'}
            ]

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
                        brand="3D MOLEKÜL GÖRÜNTÜLEME",
                        brand_href=HttpResponseRedirect(reverse("bioinformatic:single_molecule_3d_view")).url,
                        color="primary",
                        dark=True,
                        brand_external_link=True,
                        sticky='top',
                        className="shadow-lg bg-body rounded mt-1 mb-1 mr-1 ml-1",
                    ),

                    dbc.Card(
                        [
                            dbc.Row(
                                [
                                    dbc.Col(
                                        [
                                            dcc.Tabs(
                                                id='mol3d-tabs', children=[
                                                    dcc.Tab(
                                                        label='AÇIKLAMA',
                                                        children=html.Div(
                                                            className='control-tab mt-2',
                                                            children=[
                                                                html.H4(["YAPIYA İLİŞKİN BİLGİLER"], className="mt-2"),
                                                                html.P(
                                                                    f"İD : {str(file.name[:4]).upper()}"),

                                                                html.P(
                                                                    f"ADI : {structure.get('name')}"
                                                                ),

                                                                html.P(
                                                                    f"Yayınlanma Tarihi : {structure.get('release_date')}"
                                                                ),

                                                                html.P(
                                                                    f"Metod : {structure.get('structure_method')}"
                                                                ),

                                                                html.P(
                                                                    f"Referans : {structure.get('journal')}",
                                                                    className="text-align-justify",
                                                                ),

                                                                html.A(['Protein Databankta Görüntüle'],
                                                                       href=f'https://www.rcsb.org/structure/{str(file.name[:4]).lower()}',
                                                                       target="_blank"),

                                                            ]
                                                        )
                                                    ),

                                                    dcc.Tab(
                                                        label='TABLO',
                                                        children=[
                                                            html.Div(
                                                                [
                                                                    dash_table.DataTable(
                                                                        id="zooming-table",
                                                                        columns=columns,
                                                                        data=df.to_dict("records"),
                                                                        row_selectable="single",
                                                                        page_size=10,
                                                                        filter_action='native',
                                                                        filter_options={"placeholder_text": "Filtrele"},
                                                                        editable=True,
                                                                        style_cell={'textAlign': 'center'},
                                                                        style_header={
                                                                            'backgroundColor': 'white',
                                                                            'fontWeight': 'bold'
                                                                        }
                                                                    ),

                                                                ], className="mx-auto"),

                                                            html.Div([
                                                                html.Button("Su molekülünü kaldır", id="remove-water",
                                                                            className='btn btn-primary mt-2 mx-auto'),
                                                            ]),

                                                            html.P(id="water-size", children=[])
                                                        ]
                                                    ),

                                                    dcc.Tab(
                                                        label='ARAÇLAR',
                                                        value='view-options',
                                                        children=[

                                                            html.Label("Görünüm", className="fw-bolder mt-2"),

                                                            dcc.Dropdown(
                                                                id='visual-type',
                                                                options=[
                                                                    {'label': 'Çubuk', 'value': 'stick'},
                                                                    {'label': 'Şerit', 'value': 'cartoon'},
                                                                    {'label': 'Küre', 'value': 'sphere'},
                                                                ],
                                                                value='cartoon',
                                                            ),

                                                            html.Label("Renklendirme Türü", className="fw-bolder mt-2"),

                                                            dcc.Dropdown(
                                                                id='color-type',
                                                                options=[
                                                                    {'label': 'Atom', 'value': 'atom'},
                                                                    {'label': 'Bölge', 'value': 'residue'},
                                                                    {'label': 'Bölge türü', 'value': 'residue_type'},
                                                                    {'label': 'Zincir', 'value': 'chain'},
                                                                ],
                                                                value='residue', className="mt-2"
                                                            ),

                                                            html.Label("Seçme Türü", className="fw-bolder mt-2"),

                                                            dcc.Dropdown(
                                                                id='select-type',
                                                                options=[
                                                                    {'label': 'Atom', 'value': 'atom'},
                                                                    {'label': 'Bölge', 'value': 'residue'},
                                                                    {'label': 'Zincir', 'value': 'chain'},
                                                                ],
                                                                value='atom', className="mt-2"
                                                            ),

                                                            html.P(["Seçtiğiniz Bölge"], className="fw-bolder mt-2"),
                                                            html.Div(
                                                                id='select-atom-output',
                                                                className="mx-auto",
                                                                children=[]
                                                            ),
                                                        ]
                                                    ),

                                                ], className="mb-2"
                                            ),
                                        ], md=4
                                    ),

                                    dbc.Col(
                                        [
                                            dash_bio.Molecule3dViewer(
                                                id="visual_output",
                                                modelData=data,
                                                styles=styles,
                                                style={'marginRight': 'auto', 'marginLeft': 'auto'},
                                                selectionType='atom',
                                                width=500,
                                                height=600
                                            ),
                                        ], md=8,
                                    ),
                                ],
                            ),
                        ], className="shadow-lg p-3 bg-body rounded mr-1 ml-1 mb-2"
                    ),
                ],
            )

            @app.callback(
                Output("visual_output", "styles"),
                Output("visual_output", "selectionType"),
                Input("visual-type", "value"),
                Input("color-type", "value"),
                Input("select-type", "value"),
            )
            def visual_update(visualization_type, color_element, select_type):
                styles = create_mol3d_style(
                    atoms=data["atoms"], visualization_type=visualization_type, color_element=color_element
                )
                return styles, select_type

            @app.callback(
                Output("visual_output", "zoomTo"),
                Output("visual_output", "labels"),
                Input("zooming-table", "selected_rows"),
                prevent_initial_call=True,
            )
            def residue_zoom(selected_row):
                row = df.iloc[selected_row]
                row['positions'] = row['positions'].apply(lambda x: [float(x) for x in x.split(',')])

                return [
                    {
                        "sel": {"chain": row["chain"], "resi": row["residue_index"]},
                        "animationDuration": 1500,
                        "fixedPath": True,
                    },

                    [
                        {
                            "text": "Bölge ADI: {}".format(row["residue_name"].values[0]),
                            "position": {
                                "x": row["positions"].values[0][0],
                                "y": row["positions"].values[0][1],
                                "z": row["positions"].values[0][2],
                            },
                        }
                    ],

                ]

            @app.callback(
                Output("zooming-table", "data"),
                Output("water-size", "children"),
                Output("visual_output", "modelData"),
                Input("remove-water", "n_clicks"),
                prevent_initial_call=True,
            )
            def remove_water(n_clicks):

                if n_clicks:
                    data_update = parser.mol3d_data()

                    atoms = [i for i in data_update["atoms"] if not "HOH" in i.get("residue_name")]

                    water_molecule = len(data_update["atoms"]) - len(atoms)

                    data_update['atoms'] = atoms

                    data_update['bonds'] = data_update['bonds'][:-water_molecule]

                    modelData = data_update

                    df = pd.DataFrame(atoms)

                    df['positions'] = df['positions'].apply(lambda x: ', '.join(map(str, x)))

                    children = f'{water_molecule} su molekülü kaldırıldı.'

                    return df.to_dict("records"), children, modelData

            @app.callback(
                Output('select-atom-output', 'children'),
                Input('visual_output', 'selectedAtomIds'),

            )
            def show_selected_atoms(atom_ids):

                if atom_ids is None or len(atom_ids) == 0:
                    return 'Henüz bir bölge seçmediniz. Molekülün üzerine tıklayın'

                return [html.Div([
                    html.Div('Element: {}'.format(data['atoms'][atm]['elem'])),
                    html.Div('Zincir: {}'.format(data['atoms'][atm]['chain'])),
                    html.Div('Bölge: {}'.format(data['atoms'][atm]['residue_name'])),
                    html.Div('Tablo sırası: {}'.format(data['atoms'][atm]['serial'])),
                    html.Br()
                ]) for atm in atom_ids].pop()

        return HttpResponseRedirect("/laboratuvarlar/bioinformatic-laboratuvari/app/molecule-3d-viewer/")

    return render(request, 'bioinformatic/form.html', {'form': form, 'title': 'Tekli 3D MOLEKÜL GÖRÜNTÜLEME'})


def multi_molecule_view(request):
    if request.user.is_anonymous:
        from django.conf import settings
        messages.error(request, "Lütfen Giriş Yapınız")
        return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))

    if BioinformaticModel.objects.filter(user=request.user, tool="multi_molecule_view").exists():
        BioinformaticModel.objects.filter(user=request.user, tool="multi_molecule_view").delete()

    external_stylesheets = [dbc.themes.BOOTSTRAP]

    app = DjangoDash('multi-molecule-viewer', external_stylesheets=external_stylesheets,
                     title='Çoklu 3D MOLEKÜL GÖRÜNTÜLEME', add_bootstrap_links=True)

    form = MultiMoleculeViewForm(request.POST or None, request.FILES or None)

    data_path = os.path.join(BASE_DIR, "media", "laboratory", f"{request.user}\\")

    if request.method == "POST":

        if form.is_valid():

            if request.FILES['files'].size > 15 * 1024 * 1024:
                messages.error(request, "Dosya boyutu en fazla 15mb olmalıdır.")
                return HttpResponseRedirect(request.path)

            files = form.cleaned_data["files"]

            obj = BioinformaticModel.objects.create(
                user=request.user,
                tool="multi_molecule_view",
            )

            file_name = []

            for file in files:
                obj.records_files.create(file=file)
                file_name.append(str(file).upper().rsplit(".PDB")[0])

            dropdown_options = [{"label": i, "value": i} for i in file_name]

            representation_options = [
                {"label": "backbone", "value": "backbone"},
                {"label": "ball+stick", "value": "ball+stick"},
                {"label": "cartoon", "value": "cartoon"},
                {"label": "hyperball", "value": "hyperball"},
                {"label": "licorice", "value": "licorice"},
                {"label": "axes+box", "value": "axes+box"},
                {"label": "helixorient", "value": "helixorient"},
                {"label": "ribbon", "value": "ribbon"},
                {"label": "rope", "value": "rope"},
                {"label": "spacefill", "value": "spacefill"},
                {"label": "surface", "value": "surface"},
                {"label": "trace", "value": "trace"},
                {"label": "tube", "value": "tube"},
                {"label": "unitcell", "value": "unitcell"},
            ]

            camera = [
                {"label": "OTOMATİK", "value": "auto"},
                {"label": "DÜŞÜK", "value": "low"},
                {"label": "ORTA", "value": "medium"},
                {"label": "YÜKSEK", "value": "high"},
            ]

            back_color = [
                {"label": "SİYAH", "value": "black"},
                {"label": "BEYAZ", "value": "white"},
            ]

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
                                className="float-right",

                            ),
                        ],
                        brand="3D MOLEKÜL GÖRÜNTÜLEME",
                        brand_href=HttpResponseRedirect(reverse("bioinformatic:multiple_molecule_3d_view")).url,
                        color="primary",
                        dark=True,
                        brand_external_link=True,
                        sticky='top',
                        className="shadow-lg bg-body rounded mt-2 mb-1 mr-1 ml-1",
                    ),

                    dbc.Card(
                        [
                            dbc.Row(
                                [
                                    dbc.Col(
                                        [
                                            dcc.Tabs(
                                                id='mol3d-tabs', children=[
                                                    dcc.Tab(
                                                        label='AÇIKLAMA',
                                                        children=html.Div(
                                                            className='control-tab mt-2',
                                                            children=[
                                                                html.H4(["YAPIYA İLİŞKİN BİLGİLER"], className="mt-2"),
                                                                html.Label(["İD : "]),
                                                                html.Span(
                                                                    [f"{str(file).upper()}," for file in file_name],
                                                                ),

                                                            ]
                                                        )
                                                    ),

                                                    dcc.Tab(
                                                        label='GÖRÜNTÜLEME',
                                                        value='view-options',
                                                        children=[

                                                            html.Label("Molekül Seçin", className="fw-bolder mt-2"),

                                                            dcc.Dropdown(
                                                                id="molecule-select",
                                                                options=dropdown_options,
                                                                placeholder="Molekül Seçin",
                                                                value=file_name,
                                                                multi=True,
                                                            ),

                                                            html.Label("Görünüm", className="fw-bolder mt-2"),

                                                            dcc.Dropdown(
                                                                id="molecule-representation",
                                                                options=representation_options,
                                                                multi=True, value=["cartoon", "axes+box"]
                                                            ),

                                                            dcc.RadioItems(
                                                                id="nglstyle-radio",
                                                                options=[
                                                                    {'label': 'Ayrı', 'value': "True"},
                                                                    {'label': 'Bağımlı', 'value': "False"},
                                                                ],
                                                                value="False",
                                                                className="mt-2",
                                                                inline=True,
                                                                labelClassName="mr-2",
                                                                inputClassName="mr-2 ml-1"
                                                            ),

                                                            html.Label("Arka Plan Rengi", className="fw-bolder mt-1"),

                                                            dcc.Dropdown(
                                                                id="ngl-stage-color-dropdown",
                                                                value='white',
                                                                options=back_color
                                                            ),

                                                            html.Label("Görüntü Kalitesi", className="fw-bolder mt-2"),

                                                            dcc.Dropdown(
                                                                id="ngl-stage-quality-dropdown",
                                                                value='auto',
                                                                options=camera
                                                            ),

                                                            html.Label("Kamera", className="fw-bolder mt-2"),

                                                            dcc.Dropdown(
                                                                id="ngl-stage-camera-dropdown",
                                                                value='perspective',
                                                                options=[
                                                                    {"label": s.capitalize(), "value": s}
                                                                    for s in ["perspective", "orthographic"]
                                                                ]
                                                            ),

                                                            html.Button(id='save-img', n_clicks=0,
                                                                        children="Görüntüyü İndir",
                                                                        className="btn btn-primary col-12 mt-2"),
                                                        ]
                                                    ),

                                                ], className="mb-2"
                                            ),

                                        ], md=3
                                    ),

                                    dbc.Col(
                                        [
                                            dash_bio.NglMoleculeViewer(id="molecule-output", height=600, width=1050),
                                        ], md=9, className="mx-auto"
                                    ),
                                ],
                            ),
                        ], className="shadow-lg p-3 bg-body rounded mr-1 ml-1"
                    ),
                ],
            )

            @app.callback(
                Output("molecule-output", 'data'),
                Output("molecule-output", "molStyles"),
                Output("molecule-output", "stageParameters"),
                Output("molecule-output", "downloadImage"),
                Output("molecule-output", "imageParameters"),
                Input("molecule-representation", "value"),
                Input("nglstyle-radio", "value"),
                Input("molecule-select", "value"),
                Input("ngl-stage-color-dropdown", "value"),
                Input("ngl-stage-quality-dropdown", "value"),
                Input("ngl-stage-camera-dropdown", "value"),
                Input("save-img", "n_clicks"),
            )
            def return_molecule(style, sidebyside, value, color, quality, cameraType, n_clicks):

                sidebyside_bool = sidebyside == "True"

                molstyles_dict = {
                    "representations": style,
                    "chosenAtomsColor": value,
                    "chosenAtomsRadius": 1,
                    "molSpacingXaxis": 100,
                    "sideByside": sidebyside_bool
                }

                stage_params = {
                    "quality": quality,
                    "backgroundColor": color,
                    "cameraType": cameraType
                }

                imageParameters = {
                    "antialias": True,
                    "transparent": True,
                    "trim": True,
                    "defaultFilename": f"{quality}"
                }

                downloadImage = False

                if n_clicks > 0:
                    downloadImage = True

                data_list = [ngl_parser.get_data(data_path=data_path, pdb_id=pdb_id, color='red',
                                                 reset_view=True, local=True)
                             for pdb_id in value]

                return data_list, molstyles_dict, stage_params, downloadImage, imageParameters

        return HttpResponseRedirect("/laboratuvarlar/bioinformatic-laboratuvari/app/multi-molecule-viewer/")

    return render(request, 'bioinformatic/form.html', {'form': form, 'title': 'Çoklu 3D MOLEKÜL GÖRÜNTÜLEME'})

