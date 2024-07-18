import pandas as pd
from dash import html, dcc, Input, Output, dash_table, State
from django.shortcuts import *
from django_plotly_dash import DjangoDash
import dash_bootstrap_components as dbc
import plotly.express as px
import io, base64, os
import plotly.figure_factory as ff


ex_data = px.data.iris()

label = [
    'bar', 'line', 'barpolar', 'box', 'candlestick', 'carpet', 'choropleth',
    'choroplethmapbox', 'cone', 'contour', 'contourcarpet', 'densitymapbox', 'funnel',
    'distplots', 'funnelarea', 'heatmap', 'heatmapgl', 'histogram', 'histogram2d',
    'histogram2dcontour', 'icicle', 'image', 'indicator', 'isosurface','mesh3d', 'ohlc',
    'parcats', 'parcoords', 'pie', 'pointcloud','sankey', 'scatter', 'scatter3d', 'scattercarpet',
    'scattergeo', 'scattergl', 'scattermapbox','scatterpolar', 'scattersmith', 'scatterternary',
    'splom', 'streamtube', 'sunburst','surface', 'table', 'treemap', 'violin', 'volume', 'waterfall'
]

graph_type = [
    {"label": f"{i}".upper(), "value": f"{i}"} for i in label
]


def parse_contents(contents, filename):
    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    if 'csv' in filename:
        # Assume that the user uploaded a CSV file
        return pd.read_csv(
            io.StringIO(decoded.decode('utf-8')))
    elif 'xls' or 'xlsx' in filename:
        # Assume that the user uploaded an excel file
        return pd.read_excel(io.BytesIO(decoded))


def files_table(request):
    external_stylesheets = [dbc.themes.BOOTSTRAP]

    app = DjangoDash('files-stats', external_stylesheets=external_stylesheets,
                     title="Dosya ve Grafik", add_bootstrap_links=True)

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
                                                label='Açıklama',
                                                children=[
                                                    html.P(
                                                        ["Dinamik istatistiksel uygulamasına hoş geldiniz."],
                                                        className="mt-2"
                                                    ),

                                                    html.P(
                                                        [
                                                            "Excel yada csv uzantılı dosyanızı seçtikten sonra girdiğiniz veriler, istatistik ve korelasyon tablonuz ile grafiğiniz oluşacaktır."
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
                                                            'borderWidth': '1px', 'borderStyle': 'solid',
                                                            'borderRadius': '5px', 'textAlign': 'center',
                                                            'margin': '10px', 'marginLeft': 3
                                                        },
                                                    ),

                                                    html.Small(['Dosya .xlsx, .xls yada .csv uzantılı olmalıdır.'],
                                                               className="text-danger ml-1"),

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

                                                    html.Small(["** Polar grafiklerde yarıçap"], className="text-danger"),

                                                    html.P("Y Ekseni", style={'font-weight': 'bold'},
                                                               className="ml-2 text-small"),

                                                    dcc.Dropdown(
                                                        id="y-axis",
                                                        className='col-11',
                                                        style={'marginLeft': -4, 'marginTop': -3},
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

                                                    html.Small(["*Scatter grafikte kullanılabilir"], className="text-danger ml-2"),

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
                                                label='Hist',
                                                children=[

                                                    html.Label("Histogram Türü", style={'font-weight': 'bold'},
                                                               className="ml-2 text-small mt-1"),

                                                    dcc.Dropdown(
                                                        id="histnorm",
                                                        className='col-11',
                                                        options=[
                                                            {'label': 'Yüzde', 'value': 'percent'},
                                                            {'label': 'Olasılık', 'value': 'probability'},
                                                            {'label': 'Olasılık-Yoğunluk', 'value': 'probability density'},
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

                                                    dcc.Checklist(
                                                        id="text_auto",
                                                        className="ml-2 text-small mt-1",
                                                        options=[
                                                            {'label': ' Yazıyı göster', 'value': "True"},
                                                        ],
                                                        value=["True"],
                                                    ),

                                                    dcc.Checklist(
                                                        id="markers",
                                                        className="ml-2 text-small mt-1",
                                                        options=[
                                                            {'label': ' İşaret', 'value': "True"},
                                                        ],
                                                        value=["True"],
                                                    ),

                                                    dcc.Checklist(
                                                        id="cumulative",
                                                        className="ml-2 text-small mt-1",
                                                        options=[
                                                            {'label': ' Kümulatif', 'value': "True"},
                                                        ],
                                                        value=["True"],
                                                    )
                                                ]
                                            ),

                                        ],
                                    ),
                                ], md=4, className="p-3 mb-3"
                            ),

                            dbc.Col(
                                [

                                    html.Label("Tablo", style={'font-weight': 'bold'}, className="mt-1"),

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
                                                    dcc.Graph(id='datatable-upload-graph')
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
        Output('datatable-upload-container', 'data'),
        Output('datatable-upload-container', 'columns'),
        Input('datatable-upload', 'contents'),
        State('datatable-upload', 'filename'),
    )
    def update_output(contents, filename):
        if contents is None:
            return [{}], []
        df = parse_contents(contents, filename)
        return df.to_dict('records'), [
            {"name": i, "id": i, 'type': 'numeric', 'deletable': True, "renamable": True, "selectable": True} for i in
            df.columns]

    @app.callback(
        Output('datatable-upload-graph', 'figure'),
        Output('stats_out', 'children'),
        Output('sel-col', 'options'),
        Output('x-axis', 'options'),
        Output('y-axis', 'options'),
        Output('color', 'options'),
        Output('st_corr', 'children'),
        Output('st_corr', 'label'),

        Input('datatable-upload-container', 'data'),
        Input('datatable-upload-container', 'columns'),
        Input('graph-title', 'value'),
        Input('graph-type', 'value'),
        Input('sel-col', 'value'),
        Input('x-axis', 'value'),
        Input('y-axis', 'value'),
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
        prevent_initial_call=True,
    )
    def display_graph(rows, columns, title, graph_type, selected_columns, x_axs, y_axs, color,
                      trend_line, marginal_x, marginal_y, corr_method, histnorm, barmode, text_auto, markers,cumulative):

        global fig
        facet_row = None

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

        return fig, stats_desc, sel_col, x, y, renk, stats_corr, f'{corr_method.capitalize()} Korelasyon'

    return HttpResponseRedirect("/laboratuvarlar/bioinformatic-laboratuvari/app/files-stats/")
