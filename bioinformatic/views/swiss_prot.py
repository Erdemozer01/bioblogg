import os
from pathlib import Path

import plotly

from bioinformatic.models import SwissProtModel
import Bio.SwissProt
from django.views import generic

from bioinformatic.forms.file import FileReadForm
from bioinformatic.forms.url import UrlForm
from django.shortcuts import *
from Bio import SwissProt, ExPASy
from urllib.request import urlopen

BASE_DIR = Path(__file__).resolve().parent.parent
path = os.path.join(BASE_DIR, 'files\\')


def handle_uploaded_file(f):
    with open(path + f.name, 'wb+') as destination:
        for chunk in f.chunks():
            destination.write(chunk)


def swiss_prot_file(request):
    global records_p, handle
    form = FileReadForm(request.POST or None, request.FILES or None)

    if request.method == "POST":

        if form.is_valid():

            try:

                handle_uploaded_file(request.FILES['file'])

                file_path = os.path.join(BASE_DIR, "files\\{}".format(form.cleaned_data['file']))

                handle = open(file_path)

                records = SwissProt.read(handle)

                authors = [ref.authors for ref in records.references]

                title = [ref.title for ref in records.references]

                handle.close()

                os.remove(file_path)

                return render(request, "bioinformatic/swiss/result.html",
                              {"file": records, "authors": authors, "title": title, "bre": "Dosya İçerği"})

            except Bio.SwissProt.SwissProtParserError:
                return render(request, 'bioinformatic/fasta/notfound.html',
                              {"msg": "Hatalı Dosya Türü, Lütfen Swiss Prot Dosyası giriniz",
                               "url": reverse("bioinformatic:swiss_prot_file"),
                               "bre": "Hata"})

            except ValueError:
                try:

                    file_path = os.path.join(BASE_DIR, "files\\{}".format(form.cleaned_data['file']))

                    handle = open(file_path)

                    records_p = SwissProt.parse(handle)

                    if SwissProtModel.objects.exists():
                        SwissProtModel.objects.all().delete()
                        for record in records_p:
                            SwissProtModel.objects.create(
                                sequence=record.sequence,
                                organism=record.organism,
                                accessions=record.accessions[0],
                                taxonomy_id=record.taxonomy_id[0],
                                sequence_length=record.sequence_length
                            )
                        records_p.close()
                        handle.close()
                        return redirect("bioinformatic:swiss_prot_list")

                    else:
                        for record in records_p:
                            SwissProtModel.objects.create(sequence=record.sequence, organism=record.organism,
                                                          accessions=record.accessions[0],
                                                          taxonomy_id=record.taxonomy_id[0],
                                                          sequence_length=record.sequence_length).save()

                except IndexError:
                    records_p.close()
                    handle.close()
                    return redirect("bioinformatic:swiss_prot_list")

            finally:
                try:
                    file_path = os.path.join(BASE_DIR, "files\\{}".format(form.cleaned_data["file"]))
                    cls = open(file_path, 'r')
                    cls.close()
                    os.remove(file_path)
                except:
                    file_path = os.path.join(BASE_DIR, "files\\{}".format(form.cleaned_data["file"]))
                    os.remove(file_path)


        else:
            form = FileReadForm()
            return reverse("bioinformatic:swiss_prot_file")

    return render(request, "bioinformatic/swiss/file.html", {"form": form, "bre": "Swiss-Prot Dosya Okuması"})


def swiss_list_view(request):
    accessions = []
    seq_len = []

    organizm = []

    for accession in SwissProtModel.objects.all():
        accessions.append(accession.accessions)

    for seq in SwissProtModel.objects.all():
        seq_len.append(int(seq.sequence_length))

    for organizma in SwissProtModel.objects.all():
        organizm.append(organizma.organism)

    import plotly.graph_objects as go
    import plotly.express as px

    fig = go.Figure(
        data=[go.Bar(y=seq_len, x=accessions, hovertext=organizm)],
        layout_title_text="Swiss-Prot"
    )

    fig2 = px.bar(x=accessions, y=seq_len,category_orders={'organizm': organizm}, color=organizm, title="Swiss-Prot")

    fig.update_layout(xaxis_title="Erişim Numarası", yaxis_title="Sekans Uzunluğu")

    fig2.update_layout(xaxis_title="Erişim Numarası", yaxis_title="Sekans Uzunluğu")

    return render(request, "bioinformatic/swiss/table.html", {'object_list': SwissProtModel.objects.all(), "fig": fig, 'fig2':fig2})


class SwissProtDetailView(generic.DetailView):
    template_name = "bioinformatic/swiss/detail.html"
    model = SwissProtModel


def swiss_prot_url(request):
    form = UrlForm(request.POST or None)

    if request.method == "POST":
        if form.is_valid():

            url = form.cleaned_data["url"]

            handle = urlopen(url)

            swiss_file = open(path + "swiss_prot.txt", "w")

            for i in handle:
                swiss_file.writelines(str(i).replace("b'", ""))

            return redirect("bioinformatic:swiss_prot_download")

    return render(request, "bioinformatic/swiss/url.html", {"form": form, "bre": "Swiss-Prot Url Dosya Okuması"})
