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
                     title="Tablo ve Grafik", add_bootstrap_links=True, update_title="Güncelleniyor...")

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
                                                label='Açıklama',
                                                children=[
                                                    html.P(
                                                        ["Dinamik istatistik uygulamasına hoş geldiniz."],
                                                        className="mt-2"
                                                    ),

                                                    html.P(
                                                        [
                                                            "Tabloyu istediğiniz gibi düzenleyip grafik ve tabloları görüntüleyebilirsiniz."

                                                        ],

                                                    ),

                                                    html.P(
                                                        [
                                                            "Görüntülemede sorun yaşarsanız sayfayı yenileyin yada ayarlarınızı ve verilerinizi gözden geçirin"
                                                        ],

                                                    ),
                                                ]
                                            ),

                                            dbc.Tab(
                                                label='Tablo',
                                                children=[

                                                    html.Label("Tablo ayarları", style={'font-weight': 'bold'},
                                                               className="ml-2 text-small mt-2"),

                                                    dcc.Input(
                                                        id='adding-rows-name',
                                                        className="form-control ml-2 col-11",
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
                                                        className="col-12",
                                                        options=["pearson", "kendall", "spearman"],
                                                        value="pearson",
                                                        style={'marginLeft': -4, 'marginTop': -9},
                                                    ),

                                                    html.Button('Tanımlayıcı istatistik tablosunu indir', id='ex_std',
                                                                n_clicks=0,
                                                                className='btn btn-sm btn-outline-primary mt-1 ml-2 col-10'),

                                                    html.Button('Korelasyon tablosunu indir', id='ex_cor', n_clicks=0,
                                                                className='btn btn-sm btn-outline-primary mt-1 ml-2 col-10'
                                                                ),

                                                    dcc.Download(id="download-std-xlsx"),
                                                    dcc.Download(id="download-cor-xlsx"),

                                                ]
                                            ),

                                            dbc.Tab(
                                                label='Grafik',
                                                children=[

                                                    html.Label("Grafik başlığı", style={'font-weight': 'bold'},
                                                               className="ml-2 mt-1"),

                                                    dcc.Input(id="graph-title",
                                                              className="form-control ml-2 col-11",
                                                              type="text",
                                                              value="Demo Grafik Başlığı"),

                                                    html.Label("Grafik türü", style={'font-weight': 'bold'},
                                                               className="ml-2 text-small mt-1"),

                                                    dcc.Dropdown(
                                                        id="graph-type",
                                                        className='col-11',
                                                        value='heatmap',
                                                        options=graph_type,
                                                        style={'marginLeft': -4},
                                                    ),

                                                    html.Label("X Ekseni", style={'font-weight': 'bold'},
                                                               className="ml-2 text-small mt-1"),

                                                    dcc.Dropdown(
                                                        id="x-axis",
                                                        className='col-11',
                                                        style={'marginLeft': -4},
                                                    ),

                                                    html.Small(["** Polar grafiklerde yarıçap"],
                                                               className="text-danger"),

                                                    html.P("Y Ekseni", style={'font-weight': 'bold'},
                                                           className="ml-2 text-small"),

                                                    dcc.Dropdown(
                                                        id="y-axis",
                                                        className='col-11',
                                                        style={'marginLeft': -4, 'marginTop': -3},
                                                    ),

                                                    html.Label("Z Ekseni", style={'font-weight': 'bold'},
                                                               className="ml-2 text-sm"),

                                                    dcc.Dropdown(
                                                        id="z-axis",
                                                        className='col-11',
                                                        style={'marginLeft': -4},
                                                    ),

                                                    html.Label("Lejant Seçiniz", style={'font-weight': 'bold'},
                                                               className="ml-2 text-small mt-1"),

                                                    dcc.Dropdown(
                                                        id="color",
                                                        className='col-11',
                                                        style={'marginLeft': -4},
                                                    ),

                                                    html.Label("Guruplandır", style={'font-weight': 'bold'},
                                                               className="ml-2 text-small mt-1"),

                                                    dcc.Dropdown(
                                                        id="sel-col",
                                                        className='col-11',
                                                        multi=True,
                                                        style={'marginLeft': -4},
                                                    ),

                                                    html.Label("Lineer Regresyon Çizgisi",
                                                               style={'font-weight': 'bold'},
                                                               className="ml-2 text-small mt-1"),

                                                    dcc.Dropdown(
                                                        id="trendline",
                                                        className='col-11',
                                                        options=[
                                                            {'label': 'Ordinary Least Squares (OLS)', 'value': 'ols'},
                                                        ],
                                                        style={'marginLeft': -4},
                                                    ),

                                                    html.Small(["*Scatter grafikte kullanılabilir"],
                                                               className="text-danger ml-2"),

                                                ]
                                            ),

                                            dbc.Tab(
                                                label='Yardımcı',
                                                children=[

                                                    html.Label("Kenarlık X", style={'font-weight': 'bold'},
                                                               className="ml-2 text-small mt-1"),

                                                    dcc.Dropdown(
                                                        id="marginal_x",
                                                        className='col-11',
                                                        options=["histogram", "rug", "box", "violin"],
                                                        placeholder="Seçiniz",
                                                        style={'marginLeft': -4},
                                                    ),

                                                    html.Label("Kenarlık Y", style={'font-weight': 'bold'},
                                                               className="ml-2 text-small mt-1"),

                                                    dcc.Dropdown(
                                                        id="marginal_y",
                                                        className='col-11',
                                                        options=["histogram", "rug", "box", "violin"],
                                                        placeholder="Seçiniz",
                                                        style={'marginLeft': -4},
                                                    ),

                                                ]
                                            ),

                                            dbc.Tab(
                                                label='Hist.',
                                                children=[

                                                    html.Label("Histogram Türü", style={'font-weight': 'bold'},
                                                               className="ml-2 text-small mt-1"),

                                                    dcc.Dropdown(
                                                        id="histnorm",
                                                        className='col-11',
                                                        options=[
                                                            {'label': 'Yüzde', 'value': 'percent'},
                                                            {'label': 'Olasılık', 'value': 'probability'},
                                                            {'label': 'Olasılık-Yoğunluk',
                                                             'value': 'probability density'},
                                                        ],
                                                        placeholder="Seçiniz",
                                                        style={'marginLeft': -4},
                                                    ),

                                                    html.Label("Histogram Fonksiyonu", style={'font-weight': 'bold'},
                                                               className="ml-2 text-small mt-1"),

                                                    dcc.Dropdown(
                                                        id="histfunc",
                                                        className='col-11',
                                                        options=[
                                                            {'label': 'Sayı', 'value': 'count'},
                                                            {'label': 'Toplam', 'value': 'sum'},
                                                            {'label': 'Ortalama', 'value': 'avg'},
                                                            {'label': 'min', 'value': 'min'},
                                                            {'label': 'max', 'value': 'max'},
                                                        ],
                                                        placeholder="Seçiniz",
                                                        style={'marginLeft': -4},
                                                    ),

                                                    html.Label("Bar modları", style={'font-weight': 'bold'},
                                                               className="ml-2 text-small mt-1"),

                                                    dcc.Dropdown(
                                                        id="barmode",
                                                        className='col-11',
                                                        options=[
                                                            {'label': 'Yığın', 'value': 'stack'},
                                                            {'label': 'İlişkili', 'value': 'relative'},
                                                            {'label': 'Grup', 'value': 'group'},
                                                        ],
                                                        placeholder="Seçiniz",
                                                        style={'marginLeft': -4},
                                                        value="relative"
                                                    ),

                                                    dbc.Checklist(
                                                        id="text_auto",
                                                        className="ml-2 text-small mt-1",
                                                        options=[
                                                            {'label': ' Yazıyı göster', 'value': "True"},
                                                        ],
                                                        inline=True,

                                                    ),

                                                    dbc.Checklist(
                                                        id="markers",
                                                        className="ml-2 text-small mt-1",
                                                        options=[
                                                            {'label': ' İşaret', 'value': "True"},
                                                        ],
                                                        inline=True,
                                                    ),

                                                    dbc.Checklist(
                                                        id="cumulative",
                                                        className="ml-2 text-small mt-1",
                                                        options=[
                                                            {'label': ' Kümulatif', 'value': "True"},
                                                        ],
                                                        inline=True,
                                                    )
                                                ]
                                            ),
                                        ],
                                    ),

                                ], md=4, className="p-3 mb-3"
                            ),

                            dbc.Col(
                                [

                                    html.P("Tablo", className="mt-2 border-bottom text-primary pb-2"),

                                    dash_table.DataTable(

                                        id='adding-rows-table',

                                        data=[],

                                        columns=[],

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

                                ], md=8, style={'maxWidth': '63%'}, className="mx-auto mt-3 mb-3"
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
                                                label='Tanımlayıcı İstatislik',
                                                children=[

                                                ],
                                            ),

                                            dbc.Tab(
                                                id="st_corr",
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

        [
            Output('adding-rows-graph', 'figure'),
            Output('sel-col', 'options'),
            Output('x-axis', 'options'),
            Output('y-axis', 'options'),
            Output('z-axis', 'options'),
            Output('color', 'options'),
            Output('stats_out', 'children'),
            Output('st_corr', 'children'),
            Output('st_corr', 'label'),
        ],

        [
            Input('adding-rows-table', 'data'),
            Input('adding-rows-table', 'columns'),
            Input('graph-type', 'value'),
            Input('graph-title', 'value'),
            Input('sel-col', 'value'),
            Input('x-axis', 'value'),
            Input('y-axis', 'value'),
            Input('z-axis', 'value'),
            Input('color', 'value'),
            Input('corr_method', 'value'),
            Input('trendline', 'value'),
            Input('marginal_x', 'value'),
            Input('marginal_y', 'value'),
            Input('histnorm', 'value'),
            Input('barmode', 'value'),
            Input('text_auto', 'value'),
            Input('markers', 'value'),
            Input('cumulative', 'value'),
        ],

        [
            State('adding-rows-table', 'columns')
        ]
    )
    def display_output(
            data, columns, graph_type, title,
            selected_columns, x_axs, y_axs, z_axs, color,
            corr_method, trendline, marginal_x, marginal_y, histnorm, barmode,
            text_auto, markers, cumulative,
            col):

        hover_data = None
        global fig
        facet_row = None

        df = pd.DataFrame(
            data=[[row.get(c['id'], None) for c in columns] for row in data],
            columns=[i['name'] for i in columns]
        )

        sel_col = [c.get('name') for c in col]
        x = [c.get('name') for c in col]
        y = [c.get('name') for c in col]
        z = [c.get('name') for c in col]
        renk = [c.get('name') for c in col]

        stats_desc = dbc.Table.from_dataframe(
            df.describe(), striped=True, bordered=True, hover=True, index=True, size='sm', responsive=True
        )

        stats_corr = dbc.Table.from_dataframe(
            df.corr(method=corr_method, numeric_only=True), striped=True, bordered=True, hover=True, index=True,
            size='sm', responsive=True
        )

        if graph_type == 'heatmap':
            if selected_columns:
                df = df[selected_columns]
            fig = px.imshow(df, text_auto=True, aspect="auto", title=title)

        elif graph_type == 'histogram2d':
            fig = px.density_heatmap(df, x=x_axs, y=y_axs, z=z_axs, title=title, text_auto=bool(text_auto),)

        elif graph_type == 'bar':
            fig = px.bar(df, x=x_axs, y=y_axs, color=color, title=title)

        elif graph_type == 'scatter':
            if selected_columns:
                facet_row = selected_columns[0]
            fig = px.scatter(df, x=x_axs, y=y_axs, color=color, title=title, trendline_scope="overall",
                             trendline=trendline, marginal_x=marginal_x, marginal_y=marginal_y, facet_row=facet_row)

        elif graph_type == 'line':
            fig = px.line(df, x=x_axs, y=y_axs, color=color, title=title, markers=bool(markers))

        elif graph_type == 'pie':
            fig = px.pie(df, values=x_axs, names=y_axs, title=title)

        elif graph_type == 'histogram':
            fig = px.histogram(df, x=x_axs, title=title, color=color, marginal=marginal_x, text_auto=bool(text_auto),
                               histnorm=histnorm, barmode=barmode, cumulative=bool(cumulative))

        elif graph_type == 'scatterpolar':
            fig = px.scatter_polar(df, r=x_axs, theta=y_axs,
                                   color=color, size=color,
                                   color_discrete_sequence=px.colors.sequential.Plasma_r, title=title)

        elif graph_type == 'barpolar':
            fig = px.bar_polar(df, r=x_axs, theta=y_axs, color=color,
                               color_discrete_sequence=px.colors.sequential.Plasma_r, title=title)

        elif graph_type == 'scatter3d':
            fig = px.scatter_3d(df.to_dict('records'), x=x_axs, y=y_axs, z=color, title=title, color=color)

        elif graph_type == 'box':
            fig = px.box(df, x=x_axs, y=y_axs, color=color, title=title)

        elif graph_type == 'funnel':
            fig = px.funnel(df, x=x_axs, y=y_axs, color=color, title=title)

        elif graph_type == 'contour':
            if selected_columns:
                facet_row = selected_columns[0]
            fig = px.density_contour(df, x=x_axs, y=y_axs, z=z_axs, color=color, title=title, marginal_x=marginal_x, marginal_y=marginal_y, facet_row=facet_row)

        return fig, sel_col, x, y, z, renk, stats_desc, stats_corr, f'{corr_method.capitalize()} Korelasyon'

    @app.callback(
        Output('download-std-xlsx', 'data'),
        Input('ex_std', 'n_clicks'),
        State('adding-rows-table', 'data'),
        prevent_initial_call=True,
    )
    def download_des_tables(n_clicks, data):
        df = pd.DataFrame(data)
        if n_clicks > 0:
            std = dcc.send_data_frame(df.describe().to_excel, "Tanımlayıcı istatistik Tablosu.xlsx")
        return std

    @app.callback(
        Output('download-cor-xlsx', 'data'),
        Input('ex_cor', 'n_clicks'),
        Input('corr_method', 'value'),
        State('adding-rows-table', 'data'),
        prevent_initial_call=True,
    )
    def download_cor_tables(n_clicks, method, data):
        df = pd.DataFrame(data)
        if n_clicks > 0:
            corr = dcc.send_data_frame(df.corr(method=method, numeric_only=True).to_excel,
                                       f"{method.capitalize()} Tablosu.xlsx")
        return corr

    return HttpResponseRedirect("/laboratuvarlar/bioinformatic-laboratuvari/app/table-create/")