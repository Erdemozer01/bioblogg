import os
from django.shortcuts import *
from django.core.files import File
from reportlab.lib import colors
from pathlib import Path
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO
from bioinformatic.forms.diagram import GenomeDiagramForm, RestrictionModelForms, RestrictionModelFormset, \
    RestrictionModelFormFactory, RestrictionModelFormFactory2
from bioinformatic.models import RestrictionModel, DiagramModel
from django.urls import reverse_lazy
from django.contrib.auth.decorators import login_required

BASE_DIR = Path(__file__).resolve().parent.parent.parent


def handle_uploaded_file(f):
    with open(os.path.join(BASE_DIR, "bioinformatic", "files", f"{f.name}"), 'wb+') as destination:
        for chunk in f.chunks():
            destination.write(chunk)


@login_required
def genome_diagram(request):
    global file_path, input_file, diagram
    form = GenomeDiagramForm()
    restric = RestrictionModelFormFactory2()
    enzyme_obj = RestrictionModel.objects.filter(user_id=request.user.pk)
    if request.method == "POST":
        form = GenomeDiagramForm(request.POST or None, request.FILES or None)

        if form.is_valid():

            try:

                if DiagramModel.objects.filter(user=request.user.pk).exists():
                    pass
                else:
                    DiagramModel.objects.create(user=request.user)

                handle_uploaded_file(request.FILES['file'])
                file = form.cleaned_data['file']
                sigil = form.cleaned_data['diagram_shape']
                format = form.cleaned_data['diagram_format']
                fragments = form.cleaned_data['fragment']

                diagram = DiagramModel.objects.get(user=request.user)

                input_file_path = os.path.join(BASE_DIR, "bioinformatic", "files", f"{file}")

                handle = open(input_file_path)

                records = SeqIO.read(handle, "genbank")

                gd_diagram = GenomeDiagram.Diagram(records.description)
                gd_track_for_features = gd_diagram.new_track(1, name="Annotated Features")
                gd_feature_set = gd_track_for_features.new_set()

                for feature in records.features:
                    if feature.type != "gene":
                        # Exclude this feature
                        continue
                    if len(gd_feature_set) % 2 == 0:
                        color = colors.tomato
                    else:
                        color = colors.lightblue
                    gd_feature_set.add_feature(feature, color=color, label=True, sigil=sigil)

                if enzyme_obj:
                    from Bio.SeqFeature import SeqFeature, FeatureLocation
                    enzymes = []
                    for enzyme in enzyme_obj:
                        enzymes.append((enzyme.site, enzyme.enzymes, colors.red))
                    for site, name, color in enzymes:
                        index = 0
                        while True:
                            index = records.seq.find(site, start=index)
                            if index == -1:
                                break
                            feature = SeqFeature(FeatureLocation(index, index + len(site)))
                            gd_feature_set.add_feature(
                                feature,
                                color=color,
                                name=name,
                                label=True,
                                label_size=10,
                                label_color=color,
                            )
                            index += len(site)

                gd_diagram.draw(
                    format=format,
                    orientation="landscape",
                    pagesize="A4",
                    fragments=fragments,
                    start=0,
                    end=len(records),
                    circle_core=0.5
                )

                output = os.path.join(BASE_DIR, 'plasmid_circular_nice.pdf')
                image = os.path.join(BASE_DIR, 'plasmid_linear.png')

                gd_diagram.write(output, "PDF")
                gd_diagram.write(image, "png")

                if diagram.out_file:
                    diagram.out_file.delete()

                with Path(output).open('rb') as out_obj:
                    diagram.out_file = File(out_obj, name='genome_diagram.pdf')
                    diagram.save()
                    out_obj.close()

                if diagram.image:
                    diagram.image.delete()

                with Path(image).open('rb') as img_obj:
                    diagram.image = File(img_obj, name='genome_diagram.png')
                    diagram.save()
                    img_obj.close()

                handle.close()
                os.remove(input_file_path)
                os.remove(output)
                os.remove(image)

                return render(request, "bioinformatic/diagram/result.html",
                              {'object': DiagramModel.objects.filter(user=request.user)})

            except ValueError:
                handle.close()
                os.remove(input_file_path)
                msg = "Çoklu Genbank Kaydı"
                return render(request, 'bioinformatic/fasta/notfound.html',
                              {"msg": msg, 'bre': 'Hata', 'url': reverse('bioinformatic:genome_diagram')})

        else:
            form = GenomeDiagramForm(request.POST or None, request.FILES or None)

    return render(request, "bioinformatic/diagram/analiz.html",
                  {'form': form, 'restric': restric, 'obj': enzyme_obj, 'bre': "Genom Diagram Oluşturma"})


@login_required
def add_enzyme(request, pk):
    if DiagramModel.objects.filter(user=request.user).exists():
        pass
    else:
        DiagramModel.objects.create(user=request.user)
    form = RestrictionModelFormset()
    if request.method == "POST":
        form = RestrictionModelFormset(request.POST)
        if form.is_valid():
            form.instance = DiagramModel.objects.get(user=request.user)
            form.save()
            return redirect('bioinformatic:genome_diagram')
    return render(request, "bioinformatic/diagram/add_enzym.html", {'form': form, 'bre': "Enzim Ekleme"})


@login_required
def update_enzyme(request, pk):
    user = DiagramModel.objects.get(user=request.user)
    form = RestrictionModelFormFactory()
    if request.method == "POST":
        form = RestrictionModelFormFactory(request.POST)
        if form.is_valid():
            form.instance = user
            form.save()
            return redirect('bioinformatic:genome_diagram')
    return render(request, "bioinformatic/diagram/enzyme_update.html", {'form': form, 'bre': 'Enzim Düzenleme'})


@login_required
def delete_enzyme(request, pk):
    enzymes = RestrictionModel.objects.get(pk=pk)
    enzymes.delete()
    return redirect('bioinformatic:genome_diagram')
