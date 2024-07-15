import pandas as pd
from dash import html, dcc, Input, Output, dash_table, State
from django.shortcuts import *
from django_plotly_dash import DjangoDash
import dash_bootstrap_components as dbc
import plotly.express as px

ex_data = px.data.iris()

label = [
    'bar', 'barpolar', 'box', 'candlestick', 'carpet', 'choropleth',
    'choroplethmapbox', 'cone', 'contour', 'contourcarpet', 'densitymapbox', 'funnel',
    'funnelarea', 'heatmap', 'heatmapgl', 'histogram', 'histogram2d',
    'histogram2dcontour', 'icicle', 'image', 'indicator', 'isosurface',
    'mesh3d', 'ohlc', 'parcats', 'parcoords', 'pie', 'pointcloud',
    'sankey', 'scatter', 'scatter3d', 'scattercarpet', 'scattergeo', 'scattergl', 'scattermapbox',
    'scatterpolar', 'scatterpolargl', 'scattersmith', 'scatterternary', 'splom', 'streamtube', 'sunburst',
    'surface', 'table', 'treemap', 'violin', 'volume', 'waterfall'
]

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
                className="shadow-lg bg-body rounded mb-2",
            ),

            dbc.Card(
                [

                    dbc.Row(
                        [
                            dbc.Col(
                                [

                                    dbc.Tabs(
                                        [
                                            dbc.Tab(
                                                label='Tablo',
                                                children=[

                                                    html.Label("Tablo ayarları", style={'font-weight': 'bold'},
                                                               className="ml-2 text-small"),

                                                    dcc.Input(
                                                        id='adding-rows-name',
                                                        className="form-control ml-2 col-10",
                                                        placeholder='Kolon adı',
                                                    ),

                                                    html.Button('Kolon ekle', id='adding-rows-button', n_clicks=0,
                                                                className='btn btn-sm btn-outline-primary mt-1 ml-2 col-4'),

                                                    html.Button('Satır ekle', id='editing-rows-button', n_clicks=0,
                                                                className='btn btn-sm btn-outline-primary mt-1 ml-2 col-4'
                                                                ),

                                                    html.P("Korelasyon", style={'font-weight': 'bold'},
                                                           className="ml-2 text-small mt-1"),

                                                    dcc.Dropdown(
                                                        id="corr_method",
                                                        className="col-11",
                                                        options=["pearson", "kendall", "spearman"],
                                                        value="pearson",
                                                        style={'marginLeft': -4, 'marginTop': -9},
                                                    ),
                                                ]
                                            ),

                                            dbc.Tab(
                                                label='Grafik',
                                                children=[

                                                    html.Label("Grafik türü", style={'font-weight': 'bold'},
                                                               className="ml-2 text-small mt-1"),

                                                    dcc.Dropdown(
                                                        id="graph-type",
                                                        className='col-11',
                                                        value='heatmap',
                                                        options=graph_type,
                                                        style={'marginLeft': -4},
                                                    ),

                                                    html.Label("Grafik başlığı", style={'font-weight': 'bold'},
                                                               className="ml-2 mt-1"),

                                                    dbc.Input(id="graph-title",
                                                              className="form-control ml-2 col-10",
                                                              type="text",
                                                              value="Demo Grafik Başlığı"),

                                                    html.Label("Kolon seçin", style={'font-weight': 'bold'},
                                                               className="ml-2 text-small mt-1"),

                                                    dcc.Dropdown(
                                                        id="sel-col",
                                                        className='col-11',
                                                        multi=True,
                                                        style={'marginLeft': -4},
                                                    ),

                                                ]
                                            ),

                                            dbc.Tab(
                                                label='Ayar',
                                                children=[

                                                    html.Label("X Ekseni", style={'font-weight': 'bold'},
                                                               className="ml-2 text-small mt-1"),

                                                    dcc.Dropdown(
                                                        id="x-axis",
                                                        className='col-11',

                                                        style={'marginLeft': -4},
                                                    ),

                                                    html.Label("Y Ekseni", style={'font-weight': 'bold'},
                                                               className="ml-2 text-small mt-1"),

                                                    dcc.Dropdown(
                                                        id="y-axis",
                                                        className='col-11',

                                                        style={'marginLeft': -4},
                                                    ),

                                                    html.Label("Renklendirme", style={'font-weight': 'bold'},
                                                               className="ml-2 text-small mt-1"),

                                                    dcc.Dropdown(
                                                        id="color",
                                                        className='col-11',

                                                        style={'marginLeft': -4},
                                                    ),

                                                ]
                                            ),

                                        ], className='ml-1'
                                    ),

                                ], md=3, className="p-2"
                            ),

                            dbc.Col(
                                [

                                    html.Label("Tablo", style={'font-weight': 'bold'}, className="mt-1"),

                                    dash_table.DataTable(
                                        id='adding-rows-table',
                                        data=ex_data.to_dict('records'),
                                        columns=[
                                            {"name": i, 'id': i, 'type': 'numeric', 'deletable': True,
                                             "renamable": True,
                                             "selectable": True} for i in ex_data.columns
                                        ],

                                        style_table={'overflowY': 'auto', 'overflowX': 'auto'},
                                        style_cell={'textAlign': 'center'},
                                        style_data={
                                            'whiteSpace': 'normal',
                                            'height': 'auto',
                                        },

                                        row_deletable=True,
                                        editable=True,
                                        page_action='native',
                                        page_size=6,

                                    ),

                                ], md=9, style={'maxWidth': '70%', 'marginLeft': "3%"}
                            ),
                        ]
                    ),

                    dbc.Row(

                        [

                            html.Hr(),

                            dbc.Col(
                                [

                                    dbc.Tabs(

                                        [
                                            dbc.Tab(
                                                label='Grafik',
                                                children=[
                                                    dcc.Graph(id='adding-rows-graph'),
                                                ],
                                            ),

                                            dbc.Tab(
                                                id="stats_out",
                                                label='İstatislik Tablosu',
                                                children=[

                                                ],
                                            ),

                                            dbc.Tab(
                                                id="st_corr",
                                                label='Korelasyon Tablosu',
                                                children=[

                                                ],
                                            ),
                                        ]
                                    ),
                                ], md=11, className="container-fluid"
                            ),
                        ]
                    ),
                ], className="shadow-lg p-4 bg-body rounded"
            ),

        ], className="shadow-lg p-4 bg-body rounded"
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
        State('adding-rows-table', 'columns'),
    )
    def update_columns(n_clicks, value, existing_columns):
        if n_clicks > 0:
            existing_columns.append(
                {'name': value, 'id': value, 'type': 'numeric', 'renamable': True, 'deletable': True})
        return existing_columns

    @app.callback(
        Output('adding-rows-graph', 'figure'),
        Output('sel-col', 'options'),
        Output('x-axis', 'options'),
        Output('y-axis', 'options'),
        Output('color', 'options'),
        Output('stats_out', 'children'),
        Output('st_corr', 'children'),

        Input('adding-rows-table', 'data'),
        Input('adding-rows-table', 'columns'),
        Input('graph-type', 'value'),
        Input('graph-title', 'value'),
        Input('sel-col', 'value'),
        Input('x-axis', 'value'),
        Input('y-axis', 'value'),
        Input('color', 'value'),
        Input('corr_method', 'value'),

        State('adding-rows-table', 'columns'),
    )
    def display_output(data, columns, select_graph, title, selected_columns, x_axs, y_axs, color, corr_method, col):

        hover_data = None

        df = pd.DataFrame(data=[[row.get(c['id'], None) for c in columns] for row in data],
                          columns=[i['name'] for i in columns])

        sel_col = [c.get('name') for c in col]
        x = [c.get('name') for c in col]
        y = [c.get('name') for c in col]
        renk = [c.get('name') for c in col]

        stats_desc = dbc.Table.from_dataframe(
            df.describe(), striped=True, bordered=True, hover=True, index=True, size='sm', responsive=True
        )

        stats_corr = dbc.Table.from_dataframe(
            df.corr(method=corr_method), striped=True, bordered=True, hover=True, index=True,
            size='sm', responsive=True
        )

        if select_graph == 'heatmap':
            if selected_columns:
                df = df[selected_columns]

            fig = px.imshow(df, title=title, text_auto=True, aspect="auto")

        elif select_graph == 'bar':
            if selected_columns:
                hover_data = selected_columns

            fig = px.bar(df, x=x_axs, y=y_axs, color=color, title=title, hover_data=hover_data)

        return fig, sel_col, x, y, renk, stats_desc, stats_corr

    return HttpResponseRedirect("/laboratuvarlar/bioinformatic-laboratuvari/app/table-create/")
