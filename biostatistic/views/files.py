import pandas as pd
from dash import html, dcc, Input, Output, dash_table, State
from django.shortcuts import *
from django_plotly_dash import DjangoDash
import dash_bootstrap_components as dbc
import plotly.express as px
import io, base64, os
import statsmodels.api as sm
import pandas as pd
from statsmodels.formula.api import ols
from patsy import desc
from patsy.highlevel import dmatrices
from sklearn.cluster import KMeans
import plotly.graph_objs as go

ex_data = px.data.iris()

label = [
    'bar', 'line', 'barpolar', 'box', 'funnel',
    'heatmap', 'histogram', 'histogram2d', 'k-means',
    'ohlc', 'parcats', 'parcoords', 'pie', 'pointcloud', 'sankey', 'scatter', 'scatter3d',
    'scattercarpet',
    'scattergeo', 'scattergl', 'scattermapbox', 'scatterpolar', 'scattersmith', 'scatterternary',
    'splom', 'streamtube', 'sunburst', 'surface', 'table', 'treemap', 'violin', 'volume', 'waterfall'
]

graph_type = [
    {"label": f"{i}".upper(), "value": f"{i}"} for i in label
]


def parse_contents(contents, filename):
    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    if 'csv' in filename:
        return pd.read_csv(
            io.StringIO(decoded.decode('utf-8')))
    elif 'xls' or 'xlsx' in filename:
        return pd.read_excel(io.BytesIO(decoded))


def files_table(request):
    external_stylesheets = [dbc.themes.BOOTSTRAP]

    app = DjangoDash('files-stats', external_stylesheets=external_stylesheets,
                     title="Dosya ve Grafik", add_bootstrap_links=True, update_title="Güncelleniyor...")

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
                brand="Dosya ve Grafik",
                brand_href=HttpResponseRedirect(reverse("biostatistic:files_table")).url,
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
                                                label='Giriş',
                                                children=[

                                                    html.P(["Uygulamaya ilişkin açıklamalar"],
                                                           style={'font-weight': 'bold'}, className="mt-1 ml-1 mb-1"),

                                                    html.P(
                                                        ["- Dinamik istatistik uygulamasına hoş geldiniz."],
                                                        className="mt-2 mb-0",
                                                    ),

                                                    html.P(
                                                        [
                                                            "- Excel yada csv uzantılı dosyanızı seçtikten sonra girdiğiniz veriler, "
                                                            "istatistik ve korelasyon tablonuz ile grafiğiniz oluşacaktır."
                                                        ], className="mb-0"

                                                    ),

                                                    html.P(
                                                        [
                                                            "- Görüntülemede sorun yaşarsanız sayfayı yenileyin yada ayarlarınızı ve verilerinizi gözden geçirin"
                                                        ], className="mb-0"

                                                    ),

                                                    html.Div(
                                                        [
                                                            "Not : Sayfayı yenilediğinizde verileniz kaybolacaktır."
                                                        ], className="text-danger text-sm mb-0",

                                                    ),
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

                                                    html.P("Korelasyon Metodu Seçiniz", style={'font-weight': 'bold'},
                                                           className="ml-2 text-small mt-1"),

                                                    dcc.Dropdown(
                                                        id="corr_method",
                                                        className="col-11",
                                                        options=["pearson", "kendall", "spearman"],
                                                        value="pearson",
                                                        placeholder="Seçiniz",
                                                        style={'marginLeft': -4, 'marginTop': -9},
                                                    ),

                                                    html.P("Dosya Seçiniz", style={'font-weight': 'bold'},
                                                           className="ml-2 text-small mt-1"),

                                                    dcc.Upload(
                                                        id='datatable-upload',
                                                        children=html.Div([
                                                            'Dosyanızı Seçin yada Sürükleyin',

                                                        ], className="mt-2"),
                                                        style={
                                                            'width': '85%', 'height': '80px', 'lineHeight': '60px',
                                                            'borderWidth': '1px', 'borderStyle': 'dashed',
                                                            'borderRadius': '5px', 'textAlign': 'center',
                                                            'margin': '10px', 'marginLeft': 3
                                                        },
                                                    ),

                                                    html.Small(['Not: Dosya .xlsx, .xls yada .csv uzantılı olmalıdır.'],
                                                               className="text-danger ml-1", style={'marginTop': -8}),

                                                    html.Div(id="filename_out", className="ml-1")
                                                ]
                                            ),

                                            dbc.Tab(
                                                label='Eksen',
                                                children=[

                                                    html.Label("X Ekseni", style={'font-weight': 'bold'},
                                                               className="ml-2 text-small mt-1"),

                                                    dcc.Dropdown(
                                                        id="x-axis",
                                                        className='col-11',
                                                        style={'marginLeft': -4},
                                                    ),

                                                    html.Div([
                                                        html.Small(["*Polar grafiklerde yarıçap"],
                                                                   className="text-sm text-danger ml-2"),
                                                    ]),

                                                    html.Label("Y Ekseni", style={'font-weight': 'bold'},
                                                               className="ml-2 text-sm"),

                                                    dcc.Dropdown(
                                                        id="y-axis",
                                                        className='col-11',
                                                        style={'marginLeft': -4},
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
                                                label='Yard.',
                                                labelClassName="text-sm",
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

                                                    html.Label("Guruplandır", style={'font-weight': 'bold'},
                                                               className="ml-2 text-small mt-1"),

                                                    dcc.Dropdown(
                                                        id="sel-col",
                                                        className='col-11',
                                                        multi=True,
                                                        style={'marginLeft': -4},
                                                    ),

                                                    html.Label("Kümeleme Sayısı", style={'font-weight': 'bold'},
                                                               className="ml-2 text-small mt-1"),
                                                    dbc.Input(id="cluster-count", type="number", value=3, className="col-10 ml-2"),

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

                                            dbc.Tab(
                                                label='Model',
                                                children=[
                                                    html.P(["İstatistiksel Moldellemeler"],
                                                           style={'font-weight': 'bold'}, className="mt-1 ml-1 mb-1"),

                                                    dcc.Dropdown(
                                                        id="tests",
                                                        className='col-11',
                                                        options=[
                                                            {'label': 'Ordinary Least Squares (OLS)', 'value': 'ols'},
                                                            {'label': 'ANOVA (one way)', 'value': 'anova_oneway'},
                                                            {'label': 'ANOVA (two way)', 'value': 'anova_twoway'},
                                                        ],
                                                        placeholder="Seçiniz",
                                                        style={'marginLeft': -4},
                                                    ),

                                                    html.P(["Bağımlı Değişken"],
                                                           style={'font-weight': 'bold'}, className="mt-1 ml-1 mb-1"),

                                                    dcc.Dropdown(
                                                        id="dep_val",
                                                        className='col-11',
                                                        placeholder="Seçiniz",
                                                        style={'marginLeft': -4},
                                                    ),

                                                    html.Small(["**Not: Bağımlı değişken sadece numeric olabilir."],
                                                               className="ml-1 p-2 text-danger"),

                                                    html.P(["Bağımsız Değişken"],
                                                           style={'font-weight': 'bold'}, className="mt-1 ml-1 mb-1"),

                                                    dcc.Dropdown(
                                                        id="in_dep_val",
                                                        className='col-11',
                                                        placeholder="Seçiniz",
                                                        style={'marginLeft': -4},
                                                        multi=True
                                                    ),

                                                    html.Small(["**Not : En fazla 4 değişken seçebilirsiniz."],
                                                               className="ml-1 p-2 text-danger")

                                                ]
                                            ),

                                        ],
                                    ),

                                ], md=4, className="p-3 mb-3"
                            ),

                            dbc.Col(
                                [

                                    html.P("Tablo", className="mt-4 border-bottom text-primary pb-2"),

                                    dash_table.DataTable(
                                        id='datatable-upload-container',
                                        style_table={'overflowY': 'auto', 'overflowX': 'auto'},
                                        style_cell={'textAlign': 'center'},
                                        style_data={
                                            'whiteSpace': 'normal',
                                            'height': 'auto',
                                        },
                                        row_deletable=True,
                                        editable=True,
                                        page_action='native',
                                        page_size=8,
                                    ),

                                ], md=8, style={'maxWidth': '63%'}, className="mx-auto mb-3"
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
                                                    dcc.Graph(id='datatable-upload-graph')
                                                ],
                                            ),

                                            dbc.Tab(
                                                id="stats_out",
                                                label='Tanımlayıcı İstatislik',
                                                children=html.Div([
                                                    html.Div(id='stats_out-container'),
                                                    html.Button('Sonuçları İndir', id='ex_std',
                                                                n_clicks=0,
                                                                className='btn btn-sm btn-outline-primary col-2 float-end mt-2'),
                                                    dcc.Download(id="download-std-xlsx"),
                                                ])
                                            ),

                                            dbc.Tab(
                                                id="st_corr",
                                                children=html.Div([
                                                    html.Div(id='st_corr-container'),
                                                    html.Button('Sonuçları İndir', id='ex_cor', n_clicks=0,
                                                                className='btn btn-sm btn-outline-primary col-2 float-end mt-2'
                                                                ),
                                                    dcc.Download(id="download-cor-xlsx"),
                                                ])
                                            ),

                                            dbc.Tab(
                                                id="test_results",
                                                label="Test Sonuçları",
                                                children=html.Div(
                                                    [
                                                        dbc.Row(
                                                            id="test_results-container",

                                                        ),

                                                    ]
                                                )
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
        Output('datatable-upload-container', 'data'),
        Output('datatable-upload-container', 'columns'),
        Input('datatable-upload', 'contents'),
        State('datatable-upload', 'filename'),
    )
    def update_output(contents, filename):
        if contents is None:
            return None, None
        df = parse_contents(contents, filename)
        return df.to_dict('records'), [
            {"name": i, "id": i, 'deletable': True, "renamable": True, "selectable": True} for i in
            df.columns]

    @app.callback(
        Output('st_corr', 'label'),
        Input('corr_method', 'value'),
    )
    def update_label(method):
        cor_met_label = f'{method} Korelasyon'.capitalize()
        return cor_met_label

    @app.callback(
        Output('datatable-upload-graph', 'figure'),
        Output('stats_out-container', 'children'),
        Output('sel-col', 'options'),
        Output('x-axis', 'options'),
        Output('y-axis', 'options'),
        Output('z-axis', 'options'),
        Output('color', 'options'),
        Output('st_corr-container', 'children'),

        Output('dep_val', 'options'),
        Output('in_dep_val', 'options'),
        Output('filename_out', 'children'),

        Input('datatable-upload-container', 'data'),
        Input('datatable-upload-container', 'columns'),
        Input('graph-title', 'value'),
        Input('graph-type', 'value'),
        Input('sel-col', 'value'),
        Input('x-axis', 'value'),
        Input('y-axis', 'value'),
        Input('z-axis', 'value'),
        Input('color', 'value'),
        Input('trendline', 'value'),
        Input('marginal_x', 'value'),
        Input('marginal_y', 'value'),
        Input('corr_method', 'value'),
        Input('histnorm', 'value'),
        Input('barmode', 'value'),
        Input('text_auto', 'value'),
        Input('markers', 'value'),
        Input('cumulative', 'value'),
        Input('histfunc', 'value'),
        Input('cluster-count', 'value'),
        State('datatable-upload', 'filename'),

        prevent_initial_call=True,
    )
    def display_graph(rows, columns, title, graph_type, selected_columns, x_axs, y_axs, z_axs, color,
                      trend_line, marginal_x, marginal_y, corr_method, histnorm, barmode, text_auto, markers,
                      cumulative, histfunc, n_clusters, filename):

        global fig, facet_col
        facet_row = None

        file_name = html.P([f'{(filename)}'], className='fw-bold')

        df = pd.DataFrame(
            data=[[row.get(c['id'], None) for c in columns] for row in rows],
            columns=[i['name'] for i in columns]
        )

        if df.empty or len(df.columns) < 1:
            return {
                'data': [{'x': [], 'y': [], 'type': 'bar'}]
            }

        sel_col = [c for c in df.columns]
        x = [c for c in df.columns]
        y = [c for c in df.columns]
        renk = [c for c in df.columns]
        z = [c for c in df.columns]

        stats_desc = dbc.Table.from_dataframe(
            df.describe(), striped=True, bordered=True, hover=True, index=True, size='sm', responsive=True
        )

        stats_corr = dbc.Table.from_dataframe(
            df.corr(numeric_only=True, method=corr_method), striped=True, bordered=True, hover=True, index=True,
            size='sm', responsive=True
        )

        if graph_type == 'heatmap':
            if selected_columns:
                df = df[selected_columns]
            fig = px.imshow(df, text_auto=True, aspect="auto", title=title)

        elif graph_type == 'bar':
            fig = px.bar(df, x=x_axs, y=y_axs, color=color, title=title)

        elif graph_type == 'histogram2d':
            fig = px.density_heatmap(df, x=x_axs, y=y_axs, z=z_axs, title=title, text_auto=bool(text_auto),
                                     marginal_x=marginal_x, marginal_y=marginal_y)

        elif graph_type == 'scatter':
            if selected_columns:
                facet_row = selected_columns[0]
            fig = px.scatter(df, x=x_axs, y=y_axs, color=color, title=title, trendline_scope="overall",
                             trendline=trend_line, marginal_x=marginal_x, marginal_y=marginal_y, facet_row=facet_row)

        elif graph_type == 'line':
            fig = px.line(df, x=x_axs, y=y_axs, color=color, title=title, markers=bool(markers))

        elif graph_type == 'pie':
            fig = px.pie(df, values=x_axs, names=y_axs, title=title)

        elif graph_type == 'histogram':
            fig = px.histogram(df, x=x_axs, title=title, color=color, marginal=marginal_x, text_auto=bool(text_auto),
                               histnorm=histnorm, barmode=barmode, cumulative=bool(cumulative), histfunc=histfunc)

        elif graph_type == 'scatterpolar':
            fig = px.scatter_polar(df, r=x_axs, theta=y_axs,
                                   color=color, size=color,
                                   color_discrete_sequence=px.colors.sequential.Plasma_r, title=title)

        elif graph_type == 'barpolar':
            fig = px.bar_polar(df, r=x_axs, theta=y_axs, color=color,
                               color_discrete_sequence=px.colors.sequential.Plasma_r, title=title)

        elif graph_type == 'scatter3d':
            fig = px.scatter_3d(df.to_dict('records'), x=x_axs, y=y_axs, z=z_axs, title=title, color=color)

        elif graph_type == 'box':
            fig = px.box(df, x=x_axs, y=y_axs, color=color, title=title)

        elif graph_type == 'funnel':
            fig = px.funnel(df, x=x_axs, y=y_axs, color=color, title=title)

        elif graph_type == 'k-means':
            km = KMeans(n_clusters=max(n_clusters, 1))
            df = df.loc[:, [x_axs, y_axs]]
            km.fit(df.values)
            df["cluster"] = km.labels_
            centers = km.cluster_centers_
            data = [
                go.Scatter(
                    x=df.loc[df.cluster == c, x_axs],
                    y=df.loc[df.cluster == c, y_axs],
                    mode="markers",
                    marker={"size": 8},
                    name="Kümeleme {}".format(c),
                )
                for c in range(n_clusters)
            ]
            data.append(
                go.Scatter(
                    x=centers[:, 0],
                    y=centers[:, 1],
                    mode="markers",
                    marker={"color": "#000", "size": 12, "symbol": "diamond"},
                    name="Kümeleme merkezi",
                )
            )

            layout = {"xaxis": {"title": x_axs}, "yaxis": {"title": y_axs}}

            fig = go.Figure(data=data, layout=layout)

            fig.update_layout(title=title)

        return fig, stats_desc, sel_col, x, y, z, renk, stats_corr, sel_col, sel_col, file_name

    @app.callback(
        Output('download-std-xlsx', 'data'),
        Input('ex_std', 'n_clicks'),
        State('datatable-upload-container', 'data'),
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
        State('datatable-upload-container', 'data'),
        prevent_initial_call=True,
    )
    def download_cor_tables(n_clicks, method, data):
        df = pd.DataFrame(data)
        if n_clicks > 0:
            corr = dcc.send_data_frame(df.corr(method=method, numeric_only=True).to_excel,
                                       f"{method.capitalize()} Tablosu.xlsx")
        return corr

    @app.callback(
        Output('test_results-container', 'children'),

        Input('tests', 'value'),
        Input('dep_val', 'value'),
        Input('in_dep_val', 'value'),

        State('datatable-upload-container', 'data'),
        prevent_initial_call=True,
    )
    def mod_st(test, dep_val, in_dep_val, rows):

        global comment
        df = pd.DataFrame(rows)
        df = df.dropna()

        try:
            formula = f'{dep_val} ~ {in_dep_val[0]} + {in_dep_val[1]} + {in_dep_val[2]} + {in_dep_val[3]}'
        except:
            try:
                formula = f'{dep_val} ~ {in_dep_val[0]} + {in_dep_val[1]} + {in_dep_val[2]}'
            except:
                try:
                    formula = f'{dep_val} ~ {in_dep_val[0]} + {in_dep_val[1]}'
                except:
                    try:
                        formula = f'{dep_val} ~ {in_dep_val[0]}'
                    except:
                        pass

        if test == 'ols':

            y, X = dmatrices(formula, data=df, return_type='dataframe')
            model = sm.OLS(y, X).fit()
            results = model.summary()

            figure = px.scatter(model.predict(X), title=f"{dep_val} ~ {in_dep_val} tahmin grafiği ", trendline="ols")

            if model.pvalues[0] > 0.05:
                comment = f"P value değeri 0.05'den büyük çıktığı için istatistiksel olarak anlamsızdır."
            elif model.pvalues[0] < 0.05:
                comment = f"P value değeri 0.05'den küçük çıktığı için istatistiksel olarak anlamlıdır."

            value = str(results) + "\n" + "=" * 64 + "\n" + f"P value : {model.pvalues[0]}" + "\n" + comment

            children = [
                dbc.Col(
                    [

                        dbc.Textarea(
                            id="results",
                            className="mt-1",
                            value=value,
                            style={'height': 400},
                        ),

                        html.Button('İndir', id='test_res',
                                    n_clicks=0,
                                    className='btn btn-sm btn-outline-primary col-2 mt-2 float-end'
                                    ),

                        dcc.Download(id="download-tests-result-ols"),
                    ], md=6
                ),

                dbc.Col(
                    [
                        dcc.Graph(id="predict_graph", figure=figure),
                    ], md=6
                ),
            ]

            return children

        elif test == 'anova_oneway':
            model = ols(formula, data=df).fit()
            results = sm.stats.anova_lm(model, typ=1)
            results = results.reset_index().rename(columns={'index': ''})
            children = [
                dbc.Col(
                    [

                        dash_table.DataTable(
                            id="results",
                            data=results.to_dict('records'),
                            columns=[{"id": i, "name": i, } for i in results.columns],
                            style_table={'overflowY': 'auto', 'overflowX': 'auto'},
                            style_cell={'textAlign': 'center'},
                            style_data={
                                'whiteSpace': 'normal',
                                'height': 'auto',
                            },
                            row_deletable=False,
                            editable=False,
                            page_action='native',
                            page_size=8,
                        ),

                        html.Button('İNDİR', id='test_res',
                                    n_clicks=0,
                                    className='btn btn-sm btn-outline-primary col-2 mt-2 float-end'
                                    ),

                        dcc.Download(id="download-tests-result"),
                    ], md=12
                ),
            ]
            return children

        elif test == 'anova_twoway':
            model = ols(formula, data=df).fit()
            results = sm.stats.anova_lm(model, typ=2)
            results = results.reset_index().rename(columns={'index': ''})
            children = [

                dbc.Col(
                    [

                        dash_table.DataTable(
                            id="results",
                            data=results.to_dict('records'),
                            columns=[{"id": i, "name": i, } for i in results.columns],
                            style_table={'overflowY': 'auto', 'overflowX': 'auto'},
                            style_cell={'textAlign': 'center'},
                            style_data={
                                'whiteSpace': 'normal',
                                'height': 'auto',
                            },
                            row_deletable=False,
                            editable=False,
                            page_action='native',
                            page_size=8,
                        ),

                        html.Button('İNDİR', id='test_res',
                                    n_clicks=0,
                                    className='btn btn-sm btn-outline-primary col-2 mt-2 float-end'
                                    ),

                        dcc.Download(id="download-tests-result"),
                    ], md=12
                ),
            ]
            return children

    @app.callback(
        Output('download-tests-result', 'data'),
        Input('test_res', 'n_clicks'),
        Input('tests', 'value'),

        State('results', 'data'),
        prevent_initial_call=True,
    )
    def download_anova_test_results(n_clicks, test, content):
        df = pd.DataFrame(content)
        if n_clicks > 0:
            return dcc.send_data_frame(df.to_excel,
                                       f"{test.capitalize()} Tablosu.xlsx")

    @app.callback(
        Output('download-tests-result-ols', 'data'),
        Input('test_res', 'n_clicks'),
        Input('tests', 'value'),

        State('results', 'value'),
        prevent_initial_call=True,
    )
    def download_test_results(n_clicks, test, content):

        if n_clicks > 0:
            return dict(content=content, filename=f"{test.capitalize()}-sonuçları.txt")

    return HttpResponseRedirect("/laboratuvarlar/bioinformatic-laboratuvari/app/files-stats/")
