import parmed
from django.shortcuts import *
from django.contrib import messages
from bioinformatic.forms import MoleculeViewForm
from bioinformatic.models import BioinformaticModel
from django_plotly_dash import DjangoDash
from dash import dcc, html, dash_table, Input, Output
import dash_bootstrap_components as dbc
import dash_bio as dashbio
from dash_bio.utils import PdbParser, create_mol3d_style
import pandas as pd
from Bio.PDB import parse_pdb_header


def molecule_view(request):
    if request.user.is_anonymous:
        from django.conf import settings
        messages.error(request, "Lütfen Giriş Yapınız")
        return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))

    if BioinformaticModel.objects.filter(user=request.user, tool="molecule_view").exists():
        BioinformaticModel.objects.filter(user=request.user, tool="molecule_view").delete()

    external_stylesheets = [dbc.themes.BOOTSTRAP]

    app = DjangoDash('molecule-3d-viewer', external_stylesheets=external_stylesheets,
                     title='3D MOLEKÜL GÖRÜNTÜLEME', add_bootstrap_links=True)

    form = MoleculeViewForm(request.POST or None, request.FILES or None)

    if request.method == "POST":

        if form.is_valid():

            if request.FILES['file'].size > 25 * 1024 * 1024:
                messages.error(request, "Dosya boyutu en fazla 25mb olmalıdır.")
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

                return redirect('bioinformatic:3d_molecule_view')

            except parmed.exceptions.FormatNotFound:

                messages.error(request, "Hatalı Dosya Formatı")

                return redirect('bioinformatic:3d_molecule_view')

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

            app.layout = html.Div(
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
                        brand_href=HttpResponseRedirect(reverse("bioinformatic:molecule_3d_view")).url,
                        color="primary",
                        dark=True,
                        brand_external_link=True,
                        sticky='top',
                        className="shadow-lg bg-body rounded mt-3",
                    ),

                    dbc.Container(
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
                                                                    f'İD : {structure.get('idcode')}'),

                                                                html.P(
                                                                    f'ADI : {structure.get('name')}'
                                                                ),

                                                                html.P(
                                                                    f'Yayınlanma Tarihi : {structure.get('release_date')}'
                                                                ),
                                                            ]
                                                        )
                                                    ),

                                                    dcc.Tab(
                                                        label='TABLO',
                                                        children=[
                                                            html.Div([
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
                                                        label='GÖRÜNTÜLEME',
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

                                                            html.P(["Seçtiğiniz Bölge"], className="fw-bolder mt-2"),
                                                            html.Div(id='select-atom-output', className="mx-auto",
                                                                     children=[]),
                                                        ]
                                                    ),

                                                ], className="mb-2"
                                            ),

                                        ], md=4
                                    ),

                                    dbc.Col(
                                        [
                                            dashbio.Molecule3dViewer(
                                                id="visual_output",
                                                modelData=data,
                                                styles=styles,
                                                style={'marginRight': 'auto', 'marginLeft': 'auto'},
                                                selectionType='atom',
                                                width=500
                                            ),
                                        ], md=8,
                                    ),
                                ],
                            ),
                        ], fluid=True, className="shadow-lg p-3 bg-body rounded"
                    ),
                ],
            )

            @app.callback(
                Output("visual_output", "zoomTo"),
                Output("visual_output", "labels"),
                Input("zooming-table", "selected_rows"),
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
                Output("visual_output", "styles"),
                Input("visual-type", "value"),
                prevent_initial_call=True,
            )
            def visual_update(value):

                styles = create_mol3d_style(
                    atoms=data["atoms"], visualization_type=value, color_element='residue'
                )

                return styles

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
                prevent_initial_call=True,
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

    return render(request, 'bioinformatic/form.html', {'form': form, 'title': '3D MOLEKÜL GÖRÜNTÜLEME'})
