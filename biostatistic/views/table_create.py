from dash import html, dcc, Input, Output, dash_table, State
from django.shortcuts import *
from django_plotly_dash import DjangoDash
import dash_bootstrap_components as dbc

label = ['bar', 'barpolar', 'box', 'candlestick', 'carpet', 'choropleth',
         'choroplethmapbox', 'cone', 'contour', 'contourcarpet', 'densitymapbox', 'funnel',
         'funnelarea', 'heatmap', 'heatmapgl', 'histogram', 'histogram2d',
         'histogram2dcontour', 'icicle', 'image', 'indicator', 'isosurface',
         'mesh3d', 'ohlc', 'parcats', 'parcoords', 'pie', 'pointcloud',
         'sankey', 'scatter', 'scatter3d', 'scattercarpet', 'scattergeo', 'scattergl', 'scattermapbox',
         'scatterpolar', 'scatterpolargl', 'scattersmith', 'scatterternary', 'splom', 'streamtube', 'sunburst',
         'surface', 'table', 'treemap', 'violin', 'volume', 'waterfall']

graph_type = [
    {"label": f"{i}".upper(), "value": f"{i}"} for i in label
]


def create_table(request):
    external_stylesheets = [dbc.themes.BOOTSTRAP]

    app = DjangoDash('table-create', external_stylesheets=external_stylesheets,
                     title="Tablo ve Grafik", add_bootstrap_links=True)

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
                        className="float-right",

                    ),
                ],
                brand="Tablo ve Grafik",
                brand_href=HttpResponseRedirect(reverse("biostatistic:table_create")).url,
                color="primary",
                dark=True,
                brand_external_link=True,
                sticky='top',
                className="shadow-lg bg-body rounded mt-2 mb-1 mr-2 ml-2",
            ),

            dbc.Card(
                [
                    dbc.Row(
                        [
                            dbc.Col(
                                [

                                    html.Label("Tablo ayarları", style={'font-weight': 'bold'}, className="ml-3 mt-3"),

                                    dcc.Input(
                                        id='adding-rows-name',
                                        className="form-control mt-2 ml-2 col-12",
                                        placeholder='Kolon adı',
                                        value='',
                                    ),

                                    html.Button('Kolon ekle', id='adding-rows-button', n_clicks=0,
                                                className='btn btn-outline-primary mt-2 ml-2 col-5'),

                                    html.Button('Satır ekle', id='editing-rows-button', n_clicks=0,
                                                className='btn btn-outline-primary mt-2 ml-2 col-5 mr-2'
                                                ),
                                    html.Label("Grafik türü", style={'font-weight': 'bold'}, className="ml-3 mt-3"),

                                    dcc.Dropdown(
                                        id="graph-type",
                                        className="ml-2",
                                        value='heatmap',
                                        options=graph_type
                                    ),

                                ], md=3,
                            ),

                            dbc.Col(
                                [
                                    html.Label("Tablo", style={'font-weight': 'bold'}, className="mt-1"),

                                    dash_table.DataTable(

                                        id='adding-rows-table',

                                        columns=[{
                                            'name': 'Column {}'.format(i),
                                            'id': 'column-{}'.format(i),
                                            'deletable': True,
                                            'renamable': True
                                        } for i in range(1, 5)],

                                        data=[
                                            {'column-{}'.format(i): (j + (i - 1) * 5) for i in range(1, 5)}
                                            for j in range(5)
                                        ],

                                        style_table={'overflowY': 'auto', 'overflowX': 'auto'},
                                        style_cell={'textAlign': 'center'},

                                        style_data={
                                            'whiteSpace': 'normal',
                                            'height': 'auto',
                                        },

                                        filter_action="native",
                                        sort_action="native",
                                        sort_mode="multi",
                                        column_selectable="single",

                                        filter_options={"placeholder_text": "Ara"},
                                        row_deletable=True,
                                        selected_rows=[],
                                        editable=True,

                                        page_action='native',
                                        page_size=10,
                                    ),

                                    html.Small(["Eklediğiniz Kolonu silerseniz diğer eklediğinizde silinecek"],
                                               className="text-small text-danger fst-italic")

                                ], md=8, className="mt-2 mx-auto"
                            ),
                        ]
                    ),

                    html.Hr(),
                    dcc.Graph(id='adding-rows-graph', className="m-2")

                ], className="shadow-lg mr-2 ml-2 mt-1"
            ),

        ], className="shadow-lg p-3 bg-body rounded"
    )

    @app.callback(
        Output('adding-rows-table', 'data'),
        Input('editing-rows-button', 'n_clicks'),
        State('adding-rows-table', 'data'),
        State('adding-rows-table', 'columns'),
    )
    def add_row(n_clicks, rows, columns):
        if n_clicks > 0:
            rows.append({c['id']: '' for c in columns})
        return rows

    @app.callback(
        Output('adding-rows-table', 'columns'),
        Input('adding-rows-button', 'n_clicks'),
        State('adding-rows-name', 'value'),
        State('adding-rows-table', 'columns'))
    def update_columns(n_clicks, value, existing_columns):
        if n_clicks > 0:
            existing_columns.append({
                'id': value, 'name': value,
                'renamable': True, 'deletable': True
            })
        return existing_columns

    @app.callback(
        Output('adding-rows-graph', 'figure'),
        Input('adding-rows-table', 'data'),
        Input('adding-rows-table', 'columns'),
        Input('graph-type', 'value'),
    )
    def display_output(rows, columns, type):
        return {
            'data': [
                {
                    'type': f'{type}',
                    'z': [[row.get(c['id'], None) for c in columns] for row in rows],
                    'x': [c['name'] for c in columns]
                }
            ]
        }

    return HttpResponseRedirect("/laboratuvarlar/bioinformatic-laboratuvari/app/table-create/")
