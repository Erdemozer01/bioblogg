import os.path, pandas, dash_table
from django.views import generic
from django.shortcuts import *
from pathlib import Path
from django.contrib import messages
from bioinformatic.forms.writing import SelectWritingForm, FileWritingForm
from bioinformatic.models.bioinformatic import BioinformaticModel, RecordModel
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from django.shortcuts import *
from django_plotly_dash import DjangoDash
from dash import dcc, html, Input, Output, State, MATCH, Patch
import dash_bootstrap_components as dbc
import dash_ag_grid as dag

BASE_DIR = Path(__file__).resolve().parent.parent.parent


def DashWriteFile(request):
    external_stylesheets = [dbc.themes.BOOTSTRAP]

    app = DjangoDash('write-file', external_stylesheets=external_stylesheets,
                     title='Dosya yazma', update_title="Güncelleniyor...", add_bootstrap_links=True)

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
                brand="2D MOLEKÜL GÖRÜNTÜLEME",
                brand_href=HttpResponseRedirect(reverse("bioinformatic:molecule_2d_view")).url,
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
                                    dcc.Tabs(
                                        id='mol3d-tabs', children=[

                                            dcc.Tab(
                                                label='Açıklama',
                                                children=html.Div(
                                                    className='control-tab mt-2',
                                                    children=[
                                                        html.P(["Dosya oluştuma"],
                                                               className="fw-bolder mt-2"),
                                                        html.P([
                                                            "Oluşturma istediğiniz dosya türünü oluştur sekmesinden siçiniz ve taloya "
                                                            "Kolon ve sutunları düzenleyip kayıtları ekleyin oluştura tıklayın."],
                                                            className="text-primary mt-2"),

                                                    ]
                                                )
                                            ),

                                            dcc.Tab(
                                                label='Oluştur',
                                                value='view-options',
                                                children=[

                                                    html.Label("Dosya türünü seçiniz",
                                                               className="text-primary fw-bolder mt-2"),

                                                    dcc.Dropdown(
                                                        id='file-type',
                                                        options=[
                                                            {'label': 'Fasta', 'value': 'fasta'},
                                                            {'label': 'Genbank', 'value': 'gb'},
                                                        ],
                                                        value='fasta',
                                                    ),

                                                    html.Hr(className="text-danger fw-bolder mt-2"),

                                                    html.Div(id="output-rec")

                                                ]
                                            ),

                                        ],
                                    ),

                                ], md=4, lg=4,
                            ),

                            dbc.Col(
                                [

                                    dash_table.DataTable(
                                        id='adding-rows-table',
                                        columns=[
                                            {'name': 'id', 'deletable': True, 'renamable': True},
                                            {'name': 'Tanım', 'deletable': True, 'renamable': True},
                                            {'name': 'Sekans', 'deletable': True, 'renamable': True},
                                        ],

                                        editable=True,
                                        row_deletable=True
                                    ),

                                    html.Button('Kayıt ekle', id='editing-rows-button', n_clicks=0,
                                                className='btn btn-primary float-right mt-2'),

                                ], md=8, className="mx-auto"
                            ),

                        ],
                    ),
                ], className="mr-2 ml-2"
            ),
        ], className="shadow-lg p-3 bg-body rounded mr-2 ml-2"
    )

    @app.callback(
        Output('adding-rows-table', 'data'),
        Input('editing-rows-button', 'n_clicks'),
        State('adding-rows-table', 'data'),
        State('adding-rows-table', 'columns'))
    def add_row(n_clicks, rows, columns):
        if n_clicks > 0:
            list(rows).append({c['id']: '' for c in columns})
        return rows

    @app.callback(
        Output('adding-rows-table', 'columns'),
        Input('adding-rows-button', 'n_clicks'),
        State('adding-rows-name', 'value'),
        State('adding-rows-table', 'columns'))
    def update_columns(n_clicks, value, existing_columns):
        if n_clicks > 0:
            existing_columns.append(
                {'name': value, 'deletable': True, 'renamable': True},
                {'name': value, 'deletable': True, 'renamable': True},
                {'name': value, 'deletable': True, 'renamable': True},
            )
        return existing_columns

    return HttpResponseRedirect("/laboratuvarlar/bioinformatic-laboratuvari/app/write-file/")


def file_writing_format_select(request, user):
    global file_format

    if request.user.is_anonymous:
        from django.conf import settings
        messages.error(request, "Lütfen Giriş Yapınız")
        return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))

    if BioinformaticModel.objects.filter(user=request.user, tool="DOSYA YAZMA").exists():
        BioinformaticModel.objects.filter(user=request.user, tool="DOSYA YAZMA").delete()

    form = SelectWritingForm(request.POST or None)

    if request.method == "POST":
        if form.is_valid():
            file_format = form.cleaned_data['writing_file_format']

            BioinformaticModel.objects.create(
                user=request.user,
                tool='DOSYA YAZMA',
                writing_file_format=file_format,
            )

            return HttpResponseRedirect(
                reverse("bioinformatic:file_writing", kwargs={'user': request.user, 'format': file_format}))

    return render(request, 'bioinformatic/form.html', {'form': form, 'title': 'DOSYA OLUŞTURMA'})


class FileWritingListView(generic.ListView):
    template_name = "bioinformatic/writing/list.html"
    model = RecordModel
    paginate_by = 10

    def get(self, request, *args, **kwargs):
        if request.user.is_anonymous:
            from django.conf import settings
            messages.error(request, "Lütfen Giriş Yapınız")
            return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))
        return super().get(request, *args, **kwargs)

    def get_queryset(self):
        return RecordModel.objects.filter(records__user=self.request.user, records__tool="DOSYA YAZMA")

    def get_context_data(self, *, object_list=None, **kwargs):
        context = super().get_context_data(**kwargs)
        context['count'] = self.object_list.count()
        context['format'] = BioinformaticModel.objects.filter(user=self.request.user,
                                                              tool="DOSYA YAZMA").last().writing_file_format
        context['title'] = f"{context['format'].upper()} KAYITLARI"
        return context


def FileWritingView(request, format, user):
    if request.user.is_anonymous:
        from django.conf import settings
        messages.error(request, "Lütfen Giriş Yapınız")
        return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))

    form = FileWritingForm(request.POST or None)

    if request.method == "POST":

        if form.is_valid():

            sekans = form.cleaned_data['sequence'].upper()

            sekans = sekans.replace(" ", "")
            sekans = sekans.replace("\n", "")
            sekans = sekans.replace("\t", "")
            sekans = sekans.replace("\r", "")
            sekans = sekans.replace("0", "")
            sekans = sekans.replace("1", "")
            sekans = sekans.replace("2", "")
            sekans = sekans.replace("3", "")
            sekans = sekans.replace("4", "")
            sekans = sekans.replace("5", "")
            sekans = sekans.replace("6", "")
            sekans = sekans.replace("7", "")
            sekans = sekans.replace("8", "")
            sekans = sekans.replace("9", "")

            object_list = BioinformaticModel.objects.filter(user=request.user, tool="DOSYA YAZMA")

            for obj in object_list:
                obj.record_content.create(
                    molecule=form.cleaned_data['molecule'],
                    molecule_id=form.cleaned_data['molecule_id'].replace("|", ""),
                    name=form.cleaned_data['name'],
                    description=form.cleaned_data['description'],
                    db_xrefs=form.cleaned_data['db_xrefs'],
                    annotations=form.cleaned_data['annotations'],
                    source=form.cleaned_data['source'],
                    organism=form.cleaned_data['organism'],
                    keywords=form.cleaned_data['keywords'],
                    accession=form.cleaned_data['accessions'],
                    sequence=str(sekans)
                )

            return HttpResponseRedirect(
                reverse("bioinformatic:file_writing_list", kwargs={'format': format, 'user': request.user}))

        else:

            form = FileWritingForm()

    return render(
        request, 'bioinformatic/writing/form.html',
        {
            'form': form,
            'title': f'{format.upper()} DOSYASI OLUŞTURMA',
            'format': format
        }
    )


def CreateFileView(request, user, format):
    if request.user.is_anonymous:
        from django.conf import settings
        messages.error(request, "Lütfen Giriş Yapınız")
        return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))

    obj_list = RecordModel.objects.filter(
        records__user=request.user,
        records__tool='DOSYA YAZMA',
        records__writing_file_format=format
    )

    file_path = os.path.join(BASE_DIR, 'media', f'{user}_{format}_file.{format}')

    SeqIO.write(

        [
            SeqRecord(

                seq=Seq(record.sequence.replace('\n', '').replace("\r", "")).replace("\t", ""),
                id=record.molecule_id,
                name=str(record.name).replace(" ", ""),
                description=record.description,
                dbxrefs=[record.db_xrefs],
                annotations={
                    'molecule_type': [record.molecule],
                    'source': [record.source],
                    'keywords': [str(record.keywords)],
                    'organism': [record.organism],
                    'accession': [record.accession],
                    'MDAT': [record.records.created],
                },

            )

            for record in obj_list

        ], file_path, format
    )

    return HttpResponseRedirect(reverse('bioinformatic:download'))


class RecordDetailView(generic.DetailView):
    template_name = "bioinformatic/writing/detail.html"
    model = RecordModel


def RecordDeleteView(request, pk):
    RecordModel.objects.get(pk=pk).delete()
    return redirect(request.META['HTTP_REFERER'])
