import dash_bio.utils.ngl_parser as ngl_parser
import parmed.exceptions
from django.views import generic
import os
from pathlib import Path
from dash import Dash, html, dash_table, dcc
from dash.dependencies import Input, Output, State
import dash_bio as dashbio
from dash_bio.utils import PdbParser, create_mol3d_style
import pandas as pd
from bioinformatic.forms.molecule import MoleculeForm, MultipleMoleculeForm
from bioinformatic.models import MolecularModel
from django.contrib import messages
from django.shortcuts import *
from django_plotly_dash import DjangoDash
from django.core.files import File
from django.conf import settings
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import parse_pdb_header
import dash_bootstrap_components as dbc

BASE_DIR = Path(__file__).resolve().parent.parent.parent


def handle_uploaded_file(file):
    with open(os.path.join(BASE_DIR, "bioinformatic", "files", f"{file}"), 'wb+') as destination:
        for chunk in file.chunks():
            destination.write(chunk)


def molecule3dviewer(request):
    if request.user.is_anonymous:
        from django.conf import settings
        messages.error(request, "Lütfen Giriş Yapınız")
        return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))
    obj = MolecularModel()
    if MolecularModel.objects.filter(user=request.user).exists():
        MolecularModel.objects.filter(user=request.user).delete()

    form = MoleculeForm(request.POST or None, request.FILES or None)

    if request.method == "POST":
        if form.is_valid():

            try:

                file = form.cleaned_data['file']
                pdb_id = form.cleaned_data['pdb_id'].lower()
                if file is not None:
                    handle_uploaded_file(request.FILES['file'])
                    file_path = os.path.join(BASE_DIR, "bioinformatic", "files", f"{request.FILES['file']}")

                    file_size = os.stat(file_path)
                    if file_size.st_size > 2000000:
                        msg = "Dosya Boyutu 2mb fazla olmamalı"
                        url = reverse('bioinformatic:molecule_analiz')
                        os.remove(file_path)

                        return render(request, 'bioinformatic/fasta/notfound.html',
                                      {"msg": msg, 'bre': 'Hata', "url": url})
                    else:
                        obj.in_file = File(request.FILES['file'], name=file)
                        obj.file_name = file
                        obj.id_name = str(file)[:4]

                else:
                    from Bio.PDB import PDBList
                    pdbl = PDBList()
                    pdir = os.path.join(BASE_DIR, "bioinformatic", "files")
                    pdbl.retrieve_pdb_file(pdb_id, pdir=pdir, overwrite=True)
                    file_path = os.path.join(BASE_DIR, "bioinformatic", "files", f"{pdb_id}.cif")
                    file_size = os.stat(file_path)
                    if file_size.st_size > 2000000:
                        msg = "Dosya Boyutu 2mb fazla olmamalı"
                        url = reverse('bioinformatic:molecule_analiz')
                        os.remove(file_path)
                        return render(request, 'bioinformatic/fasta/notfound.html',
                                      {"msg": msg, 'bre': 'Hata', "url": url})
                    else:
                        obj.in_file = File(Path(file_path).open('r'), name=f"{pdb_id}.cif")
                        obj.file_name = f"{pdb_id}.cif"
                        obj.id_name = str(pdb_id)[:4]

            except FileNotFoundError:
                msg = "Dosya Bulunamadı"
                url = reverse('bioinformatic:molecule_analiz')
                return render(request, 'bioinformatic/fasta/notfound.html',
                              {"msg": msg, 'bre': 'Hata', "url": url})

            obj.user = request.user
            obj.save()
            os.remove(file_path)
            object = MolecularModel.objects.filter(user=request.user).latest('created')
            return HttpResponseRedirect(
                reverse('bioinformatic:molecule3dviewer',
                        args=(request.user, object.pk, obj.id_name, obj.created.date())))

    return render(request, "bioinformatic/molecule/input.html", {'form': form, 'bre': "3D Tekli Molekül Görüntüleme"})


class Molecule3DView(generic.DetailView):
    template_name = "bioinformatic/molecule/molecule3dview.html"
    model = MolecularModel

    def get(self, request, *args, **kwargs):

        try:
            if request.user.is_anonymous:
                from django.conf import settings
                messages.error(request, "Lütfen Giriş Yapınız")
                return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))
            obj = MolecularModel.objects.filter(user=self.request.user).latest('created')
            app = DjangoDash(
                name='Molecule3dViewer',
            )

            try:
                parser = PdbParser(obj.in_file.path)
                if obj.file_name.endswith('.pdb'):
                    pbd_handle = PDBParser(PERMISSIVE=1)

                    structure = pbd_handle.get_structure(obj.id_name, obj.in_file.path)

                    obj.name = structure.header["name"]
                    obj.author = structure.header["author"]
                    obj.id_code = structure.header["idcode"]
                    obj.keywords = structure.header["keywords"]
                    obj.head = structure.header["head"]
                    obj.save()

                elif obj.file_name.endswith('.cif'):
                    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
                    mmcif_dict = MMCIF2Dict(obj.in_file.path)

                    try:
                        obj.name = mmcif_dict["_entity_src_gen.pdbx_gene_src_scientific_name"][0]
                        obj.author = mmcif_dict["_audit_author.name"][0:30]

                        if mmcif_dict["data_"]:
                            obj.id_code = mmcif_dict["data_"]
                        else:
                            obj.id_code = mmcif_dict["_pdbx_entry_details.entry_id"]
                        obj.keywords = mmcif_dict["_struct_keywords.text"][0]
                        obj.head = mmcif_dict["_struct.title"][0]
                        obj.save()
                    except:
                        pass

            except ValueError:
                messages.error(self.request, "Sayfayı yenilediğiniz İçin veriler kaybolmuştur")
                return redirect('bioinformatic:molecule_analiz')

            except parmed.exceptions.FormatNotFound:
                messages.error(self.request, "Hatalı Dosya Formatı")
                return redirect('bioinformatic:molecule_analiz')

            data = parser.mol3d_data()

            styles = create_mol3d_style(
                data['atoms'], visualization_type='cartoon', color_element='residue'
            )

            df = pd.DataFrame(data["atoms"])
            df = df.drop_duplicates(subset=['residue_name'])
            df['positions'] = df['positions'].apply(lambda x: ', '.join(map(str, x)))

            app.layout = html.Div(
                [

                    dash_table.DataTable(
                        id="zooming-specific-residue-table",
                        columns=[{"name": i, "id": i} for i in df.columns],
                        data=df.to_dict("records"),
                        row_selectable="single",
                        page_size=10,
                    ),

                    dashbio.Molecule3dViewer(
                        id="zooming-specific-molecule3d-zoomto",
                        modelData=data,
                        styles=styles,
                        style={'marginRight': 'auto', 'marginLeft': 'auto'},
                        width=950, height=600
                    ),
                ], className="container"
            )

            @app.callback(
                Output("zooming-specific-molecule3d-zoomto", "zoomTo"),
                Output("zooming-specific-molecule3d-zoomto", "labels"),
                Input("zooming-specific-residue-table", "selected_rows"),
                prevent_initial_call=True
            )
            def residue(selected_row):
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
                            "text": "Bölge Adı: {}".format(row["residue_name"].values[0]),
                            "position": {
                                "x": row["positions"].values[0][0],
                                "y": row["positions"].values[0][1],
                                "z": row["positions"].values[0][2],
                            },
                        }
                    ],
                ]

        except MolecularModel.DoesNotExist:
            redirect('bioinformatic:molecule_analiz')

        return super(Molecule3DView, self).get(request, *args, **kwargs)

    def get_context_data(self, **kwargs):
        context = super(Molecule3DView, self).get_context_data(**kwargs)
        context['bre'] = "Molekül Görüntüleme"
        context['object'] = MolecularModel.objects.filter(user=self.request.user).latest('created')
        return context


def multiplemoleculeviewer(request):
    if request.user.is_anonymous:
        from django.conf import settings
        messages.error(request, "Lütfen Giriş Yapınız")
        return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))

    if MolecularModel.objects.filter(user=request.user).exists():
        MolecularModel.objects.filter(user=request.user).delete()

    form = MultipleMoleculeForm(request.POST or None, request.FILES or None)

    if request.user.is_anonymous:
        from django.conf import settings
        messages.error(request, "Lütfen Giriş Yapınız")
        return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))

    if request.method == "POST":
        if form.is_valid():
            files = request.FILES.getlist('in_file')

            for file in files:
                MolecularModel.objects.create(
                    in_file=file, user=request.user, id_name=str(file)[:4].lower()
                )

            object = MolecularModel.objects.filter(user=request.user).latest('created')
            return HttpResponseRedirect(
                reverse('bioinformatic:ngl_molecule_detail',
                        args=(request.user, object.pk, object.created.date())))

    return render(request, "bioinformatic/molecule/input.html", {'form': form, 'bre': "Çoklu Molekül Görüntüleme"})


class MultipleMoleculeView(generic.DetailView):
    template_name = "bioinformatic/molecule/multiple.html"
    model = MolecularModel

    def get(self, request, *args, **kwargs):
        try:
            if self.request.user.is_anonymous:
                from django.conf import settings
                messages.error(request, "Lütfen Giriş Yapınız")
                return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))

            obj = MolecularModel.objects.filter(user=self.request.user)

            app = DjangoDash(
                name='MultipleMoleculeView',
            )

            pdb_id = []

            data_path = os.path.join(BASE_DIR, "media", "molecule", f"{self.request.user}\\").replace('\\', '/')

            for file in obj:
                pdb_id.append(file.id_name)

            representation_options = [
                {"label": "Omurga", "value": "backbone"},
                {"label": "Molekül Yapısı", "value": "ball+stick"},
                {"label": "Katı", "value": "cartoon"},
                {"label": "Hiperbol", "value": "hyperball"},
                {"label": "licorice", "value": "licorice"},
                {"label": "Çerçeve", "value": "axes+box"},
                {"label": "helixorient", "value": "helixorient"}
            ]

            app.layout = html.Div([

                html.H3('Gösterim Şekli: '),

                dcc.Dropdown(id="nglstyle-dropdown", options=representation_options,
                             multi=True, value=["cartoon", "axes+box"]),

                html.H3('Pozisyon:'),

                dcc.RadioItems(
                    id="nglstyle-radio",
                    options=[
                        {'label': 'Yan Yana', 'value': "True"},
                        {'label': 'Bağımsız', 'value': "False"},
                    ],
                    value="False"
                ),

                html.H3('Görüntü Boyutunu Ayarlama:'),
                html.P('Yükseklik:'),

                dcc.Slider(
                    id='height-ngl-h',
                    min=300,
                    max=1950,
                    value=600,
                    step=100,
                    marks={300: '300px', 1950: '1950px'}
                ),

                html.P('Genişlik:'),

                dcc.Slider(
                    id='width-ngl-w',
                    min=300,
                    max=1950,
                    value=600,
                    step=100,
                    marks={300: '300px', 1950: '1950px'}
                ),

                html.P("Bir dosya adı girerek görüntüyü indirebilirsiniz"),

                html.Button(id='save-ngl-simg', n_clicks=0, children="Kaydet", className="btn"),
                dcc.Input(id='file-ngl-simg', placeholder="Dosya Adı"),

                dashbio.NglMoleculeViewer(id="nglstyle-ngl"),
            ], className="container")

            @app.callback(
                Output("nglstyle-ngl", 'data'),
                Output("nglstyle-ngl", "molStyles"),
                Output("nglstyle-ngl", "height"),
                Output("nglstyle-ngl", "width"),
                Output("nglstyle-ngl", "downloadImage"),
                Output("nglstyle-ngl", "imageParameters"),
                Input("nglstyle-dropdown", "value"),
                Input("nglstyle-radio", "value"),
                Input("height-ngl-h", "value"),
                Input("width-ngl-w", "value"),
                Input("save-ngl-simg", "n_clicks"),
                State("file-ngl-simg", "value")
            )
            def return_molecule(style, sidebyside, height, width, n_clicks, filename):

                imageParameters = {
                    "antialias": True,
                    "transparent": True,
                    "trim": True,
                    "defaultFilename": filename
                }
                downloadImage = False

                if n_clicks > 0:
                    downloadImage = True

                sidebyside_bool = sidebyside == "True"

                molstyles_dict = {
                    "representations": style,
                    "chosenAtomsColor": "red",
                    "chosenAtomsRadius": 1,
                    "molSpacingXaxis": 100,
                    "sideByside": sidebyside_bool
                }

                data_list = [
                    ngl_parser.get_data(
                        data_path=data_path,
                        pdb_id=molecule,
                        color='red',
                        reset_view=True,
                        local=True
                    )
                    for molecule in pdb_id
                ]

                return data_list, molstyles_dict, height, width, downloadImage, imageParameters

        except MolecularModel.DoesNotExist:
            redirect('bioinformatic:molecule_analiz')

        return super(MultipleMoleculeView, self).get(request, *args, **kwargs)

    def get_context_data(self, **kwargs):
        context = super(MultipleMoleculeView, self).get_context_data(**kwargs)
        context['bre'] = "Çoklu Molekül Görüntüleme"
        return context


class MoleculeDetailView(generic.DetailView):
    template_name = "bioinformatic/molecule/multiple.html"
    model = MolecularModel

    def get(self, request, *args, **kwargs):

        try:

            if request.user.is_anonymous:
                from django.conf import settings
                messages.error(request, "Lütfen Giriş Yapınız")
                return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))

            obj = MolecularModel.objects.filter(user=self.request.user)
            app = DjangoDash(
                name='MultipleMoleculeView',
            )

            pdb_id = []

            data_path = os.path.join(BASE_DIR, "media", "molecule", f"{self.request.user}\\").replace('\\', '/')

            for file in obj:
                pdb_id.append(file.id_name)

            representation_options = [
                {"label": "Omurga", "value": "backbone"},
                {"label": "Molekül Yapısı", "value": "ball+stick"},
                {"label": "Katı", "value": "cartoon"},
                {"label": "Hiperbol", "value": "hyperball"},
                {"label": "licorice", "value": "licorice"},
                {"label": "Çerçeve", "value": "axes+box"},
                {"label": "helixorient", "value": "helixorient"}
            ]

            app.layout = html.Div(
                [
                    html.H3('Gösterim Şekli: '),

                    dcc.Dropdown(id="nglstyle-dropdown", options=representation_options,
                                 multi=True, value=["cartoon", "axes+box"]),

                    html.H3('Pozisyon:'),

                    dcc.RadioItems(
                        id="nglstyle-radio",
                        options=[
                            {'label': 'Yan Yana', 'value': "True"},
                            {'label': 'Bağımsız', 'value': "False"},
                        ],
                        value="False"
                    ),

                    html.H3('Görüntü Boyutunu Ayarlama:'),
                    html.P('Yükseklik:'),

                    dcc.Slider(
                        id='height-ngl-h',
                        min=300,
                        max=1950,
                        value=600,
                        step=100,
                        marks={300: '300px', 1950: '1950px'}
                    ),

                    html.P('Genişlik:'),

                    dcc.Slider(
                        id='width-ngl-w',
                        min=300,
                        max=1950,
                        value=600,
                        step=100,
                        marks={300: '300px', 1950: '1950px'}
                    ),

                    html.P("Bir dosya adı girerek görüntüyü indirebilirsiniz"),

                    html.Button(id='save-ngl-simg', n_clicks=0, children="Kaydet"),
                    dcc.Input(id='file-ngl-simg', placeholder="Dosya Adı"),

                    dashbio.NglMoleculeViewer(id="nglstyle-ngl"),

                ]
            )

            @app.callback(
                Output("nglstyle-ngl", 'data'),
                Output("nglstyle-ngl", "molStyles"),
                Output("nglstyle-ngl", "height"),
                Output("nglstyle-ngl", "width"),
                Output("nglstyle-ngl", "downloadImage"),
                Output("nglstyle-ngl", "imageParameters"),
                Input("nglstyle-dropdown", "value"),
                Input("nglstyle-radio", "value"),
                Input("height-ngl-h", "value"),
                Input("width-ngl-w", "value"),
                Input("save-ngl-simg", "n_clicks"),
                State("file-ngl-simg", "value")
            )
            def return_molecule(style, sidebyside, height, width, n_clicks, filename):

                imageParameters = {
                    "antialias": True,
                    "transparent": True,
                    "trim": True,
                    "defaultFilename": filename
                }
                downloadImage = False

                if n_clicks > 0:
                    downloadImage = True

                sidebyside_bool = sidebyside == "True"

                molstyles_dict = {
                    "representations": style,
                    "chosenAtomsColor": "red",
                    "chosenAtomsRadius": 1,
                    "molSpacingXaxis": 100,
                    "sideByside": sidebyside_bool
                }

                data_list = [
                    ngl_parser.get_data(
                        data_path=data_path,
                        pdb_id=molecule,
                        color='red',
                        reset_view=True,
                        local=True
                    )
                    for molecule in pdb_id
                ]

                return data_list, molstyles_dict, height, width, downloadImage, imageParameters

        except MolecularModel.DoesNotExist:
            redirect('bioinformatic:molecule_analiz')

        return super(MoleculeDetailView, self).get(request, *args, **kwargs)

    def get_context_data(self, **kwargs):
        context = super(MoleculeDetailView, self).get_context_data(**kwargs)
        context['bre'] = "Molekül Ayrıntısı"
        return context
