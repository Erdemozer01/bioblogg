import Bio
from django.http import HttpResponseRedirect
from django.shortcuts import render, redirect, reverse
from bioinformatic.forms.writing import FastaWritingForm
from bioinformatic.forms.add import AddFastaData
from Bio import SeqIO
from pathlib import Path
import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from bioinformatic.forms.translation import DNAFastaFileTranslateForm
from bioinformatic.forms.file import MultipleUploadFileForm, FileReadForm
from Bio.SeqUtils import GC
import pylab
from bioinformatic.models import GraphicModels, FastaRead
from django.core.files import File
from django.contrib.auth.decorators import login_required
from django.views import generic
from django.contrib import messages
import pandas as pd
import plotly.express as px
from django_plotly_dash import DjangoDash
from dash import html, dcc

BASE_DIR = Path(__file__).resolve().parent.parent
path = os.path.join(BASE_DIR, 'files\\')


def handle_uploaded_file(f):
    with open(os.path.join(BASE_DIR, "files", f"{f}"), 'wb+') as destination:
        for chunk in f.chunks():
            destination.write(chunk)


def fasta_reading(request):
    global read

    if request.user.is_anonymous:
        from django.conf import settings
        messages.error(request, "Lütfen Giriş Yapınız")
        return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))

    if FastaRead.objects.filter(user=request.user).exists():
        FastaRead.objects.filter(user=request.user).delete()

    form = FileReadForm(request.POST or None, request.FILES or None)

    if request.method == "POST":
        if form.is_valid():
            handle_uploaded_file(request.FILES['file'])
            org_type = form.cleaned_data['org_type']
            file_path = os.path.join(BASE_DIR, "files", f"{request.FILES['file']}")
            handle = open(file_path)
            try:
                if os.path.getsize(file_path) <= 5242880:
                    read = SeqIO.parse(handle, 'fasta')
                    for i in read:
                        FastaRead.objects.create(
                            user=request.user, name=i.name, sequence=i.seq,
                            protein=i.seq.translate(table=org_type).replace("*", ""), gc=GC(i.seq).__round__(2)
                        )
                else:
                    handle.close()
                    os.remove(file_path)
                    messages.error(request, "Dosya boyutu 5mb dan büyüktür.")
                    return redirect(request.META['HTTP_REFERER'])
            finally:
                handle.close()
                os.remove(file_path)

            obj = FastaRead.objects.filter(user=request.user)

            data = pd.DataFrame({
                'name': [i.name for i in obj],
                'dna_seq': [i.sequence for i in obj],
                'pro_seq': [i.protein for i in obj],
                'dna_seq_len': [len(i.sequence) for i in obj],
                'pro_seq_len': [len(i.protein) for i in obj],
                'gc': [i.gc for i in obj],
            })

            df = pd.DataFrame(

                {
                    'DNA': data["dna_seq_len"],
                    'PROTEİN': data['pro_seq_len'],
                    '%GC': data['gc']
                },

            )

            scatter_map = DjangoDash('linear')

            scatter_map.layout = html.Div([
                dcc.Graph(figure=px.scatter(
                    data_frame=data, x="dna_seq_len", y="pro_seq_len", color="name",
                    title="DNA SEKANS - PROTEİN SEKANS UZUNLUKLARI", trendline="ols",
                    labels={"dna_seq_len": "DNA SEKANS UZUNLUĞU",
                            "pro_seq_len": "PROTEİN SEKANS UZUNLUĞU",
                            'name': 'Canlı'
                            },
                )),
            ])

            return render(request, "bioinformatic/fasta/result.html",
                          {
                              'table_cov': df.cov().to_html(
                                  classes="table shadow-soft rounded text-center",
                                  border=False, header="Kovaryans", justify="center"
                              ),

                              'table_cor': df.corr().to_html(
                                  classes="table shadow-soft rounded text-center",
                                  border=False, header="Korelasyon", justify="center"
                              ),

                              'object_list': obj,
                              'org_count': obj.count(),
                              'sum_seq': sum(data['dna_seq_len']),
                              'sum_pro': sum(data['pro_seq_len']),
                              'sum_gc': sum(data['gc']),
                              'mean_pro': round(data['pro_seq_len'].mean(), 2),
                              'mean_seq': round(data['dna_seq_len'].mean(), 2),
                              'mean_gc': round(data['gc'].mean(), 2),
                              'median_pro': round(data['pro_seq_len'].median(), 2),
                              'median_seq': round(data['dna_seq_len'].median(), 2),
                              'median_gc': round(data['gc'].median(), 2),
                              'var_pro': round(data['pro_seq_len'].var(), 2),
                              'var_seq': round(data['dna_seq_len'].var(), 2),
                              'var_gc': round(data['gc'].var(), 2),
                              'std_pro': round(data['pro_seq_len'].std(), 2),
                              'std_seq': round(data['dna_seq_len'].std(), 2),
                              'std_gc': round(data['gc'].std(), 2),
                          })

    return render(request, "bioinformatic/fasta/read.html", {'form': form})


class FastaReadingDetailView(generic.DetailView):
    template_name = "bioinformatic/fasta/detail.html"
    model = FastaRead

    def get_queryset(self):
        return FastaRead.objects.filter(user=self.request.user)


@login_required
def fasta_gc_plot(request):
    form = FileReadForm(request.POST or None, request.FILES or None)
    if request.method == "POST":
        if form.is_valid():
            handle_uploaded_file(request.FILES['file'])
            file_path = os.path.join(BASE_DIR, "files", f"{request.FILES['file']}")
            image_path = os.path.join(BASE_DIR, "files", "{}_fasta_gc_plot.png".format(request.user))
            handle = open(file_path)
            read = SeqIO.parse(handle, 'fasta')
            gc_values = sorted(GC(rec.seq) for rec in read)
            pylab.plot(gc_values)
            pylab.title(
                "%GC Plot"
            )
            pylab.xlabel("Gen Sayısı")
            pylab.ylabel("%GC")
            pylab.savefig(image_path)

            obj = GraphicModels()
            if GraphicModels.objects.filter(user=request.user, graph_type="GC-PLOT").exists():
                GraphicModels.objects.filter(user=request.user, graph_type="GC-PLOT").delete()

            with Path(image_path).open('rb') as image_obj:
                obj.user = request.user
                obj.graph_type = "GC-PLOT"
                obj.fasta_gc_plot = File(image_obj, name="fasta_gc_plot.png")
                obj.save()

            handle.close()
            os.remove(file_path)
            image_handle = open(image_path)
            image_handle.close()
            os.remove(image_path)

            return HttpResponseRedirect(
                reverse('bioinformatic:fasta_gc_plot_result',
                        args=(obj.user, obj.pk, obj.graph_type.lower(), obj.created.date())))

    return render(request, "bioinformatic/fasta/read.html", {'form': form})


class GCPlotDetailView(generic.DetailView):
    template_name = "bioinformatic/fasta/gc_plot_result.html"
    model = GraphicModels

    def get_context_data(self, **kwargs):
        context = super(GCPlotDetailView, self).get_context_data(**kwargs)
        context['bre'] = "%GC Plot Sonuçları"
        return context


@login_required
def fasta_dot_plot(request):
    form = FileReadForm(request.POST or None, request.FILES or None)
    if request.method == "POST":
        if form.is_valid():
            try:
                handle_uploaded_file(request.FILES['file'])
                file_path = os.path.join(BASE_DIR, "files", f"{request.FILES['file']}")
                image_path = os.path.join(BASE_DIR, "files", "{}_fasta_dot_plot.png".format(request.user))
                handle = open(file_path)
                record = SeqIO.parse(handle, 'fasta')

                rec_one = next(record)
                rec_two = next(record)
                window = 7
                seq_one = rec_one.seq.upper()
                seq_two = rec_two.seq.upper()
                data = [
                    [
                        (seq_one[i: i + window] != seq_two[j: j + window])
                        for j in range(len(seq_one) - window)
                    ]
                    for i in range(len(seq_two) - window)
                ]

                pylab.gray()
                pylab.imshow(data)
                pylab.xlabel("%s (length %i bp)" % (rec_one.id, len(rec_one)))
                pylab.ylabel("%s (length %i bp)" % (rec_two.id, len(rec_two)))
                pylab.title("Dot plot")
                pylab.savefig(image_path)
                obj = GraphicModels()
                if GraphicModels.objects.filter(user=request.user, graph_type="Dot Plot").exists():
                    GraphicModels.objects.filter(user=request.user, graph_type="Dot Plot").delete()

                with Path(image_path).open('rb') as file_obj:
                    obj.user = request.user
                    obj.graph_type = "Dot Plot"
                    obj.fasta_dot_plot = File(file_obj, name="fasta_dot_plot.png")
                    obj.save()

                handle.close()
                os.remove(file_path)
                image_handle = open(image_path)
                image_handle.close()
                os.remove(image_path)

                return HttpResponseRedirect(
                    reverse('bioinformatic:fasta_dot_plot_result',
                            args=(obj.user, obj.pk, obj.graph_type.lower(), obj.created.date())))

            except RuntimeError:
                return redirect('bioinformatic:fasta_dot_plot')

    return render(request, "bioinformatic/fasta/read.html", {'form': form, "bre": "Fasta Dot Plot"})


class DotPlotDetailView(generic.DetailView):
    template_name = "bioinformatic/fasta/dot_plot_result.html"
    model = GraphicModels

    def get_context_data(self, **kwargs):
        context = super(DotPlotDetailView, self).get_context_data(**kwargs)
        context['bre'] = "Dot Plot Sonuçları"
        return context


def fasta_writing(request):
    form = FastaWritingForm(request.POST or None)

    if request.user.is_anonymous:
        from django.conf import settings
        messages.error(request, "Lütfen Giriş Yapınız")
        return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))

    if request.method == "POST":

        if form.is_valid():

            fasta_id = form.cleaned_data["id"]

            name = form.cleaned_data["name"]

            description = form.cleaned_data["description"]

            sequence = form.cleaned_data["sequence"].replace("\n", "").replace("\r", "")

            file_path = os.path.join(BASE_DIR, "files", f"{request.user.username}.fasta")

            record = SeqRecord(
                seq=Seq(sequence),
                id=fasta_id.encode().decode(encoding="utf-8", errors="ignore"),
                description=description,
                name=name
            )

            SeqIO.write(record, file_path, 'fasta')

            return HttpResponseRedirect(reverse('bioinformatic:fasta_download'))

        else:

            msg = "Bir hata meydana geldi"

            return render(request, 'exception/page-404.html', {
                "msg": msg
            })

    return render(request, "bioinformatic/fasta/read.html", {
        "form": form,
    })


def append_new_line(file_name, text_to_append):
    """Append given text as a new line at the end of file"""
    # Open the file in append & read mode ('a+')
    with open(file_name, "a+") as file_object:
        # Move read cursor to the start of file.
        file_object.seek(0)
        # If file is not empty then append '\n'
        data = file_object.read(100000)
        if len(data) > 0:
            file_object.write("\n")
        # Append text at the end of file
        file_object.write(text_to_append)


def fasta_add(request):
    global file_fasta
    form = AddFastaData(request.POST or None, request.FILES or None)
    if request.method == "POST":

        if form.is_valid():

            try:
                handle_uploaded_file(request.FILES["file"])
                input_file = form.cleaned_data['file']
                fasta_id = form.cleaned_data["fasta_id"]
                description = form.cleaned_data['description']
                sequence = form.cleaned_data["sequence"].replace("\n", "").replace("\r", "")
                file_fasta = os.path.join(BASE_DIR, "files", f"{input_file}")

                record = SeqRecord(
                    seq=Seq(sequence),
                    id=fasta_id.encode().decode(encoding="utf-8", errors="ignore"),
                    description=description
                ).format("fasta")

                read_input = open(file_fasta, "r").read()

                if "LOCUS" in read_input:
                    os.remove(file_fasta)
                    msg = "Lütfen Fasta Dosyası Seçiniz"
                    return render(request, "bioinformatic/fasta/notfound.html", {
                        "msg": msg
                    })

                if "#NEXUS" in read_input:
                    os.remove(file_fasta)
                    msg = "Lütfen Fasta Dosyası Seçiniz."
                    return render(request, "bioinformatic/fasta/notfound.html", {
                        "msg": msg
                    })

                append_new_line(file_fasta, str(record))
                os.rename(file_fasta, os.path.join(BASE_DIR, "files", f"{request.user}.fasta"))
                return redirect("bioinformatic:fasta_download")

            except UnicodeDecodeError:
                os.remove(file_fasta)
                msg = "Lütfen Fasta Dosyası Seçiniz."
                return render(request, "bioinformatic/fasta/notfound.html", {
                    "msg": msg
                })

        else:

            msg = "Bir hata meydana geldi"

            return render(request, "bioinformatic/fasta/notfound.html", {
                "msg": msg
            })

    return render(request, "bioinformatic/fasta/add.html", {
        "form": form,
        "bre": "Fasta Dosyası Veri Ekleme"
    })


def fasta_file_translate(request):
    global record, input_fasta_path
    form = DNAFastaFileTranslateForm(request.POST or None, request.FILES or None)
    if request.user.is_anonymous:
        from django.conf import settings
        messages.error(request, "Lütfen Giriş Yapınız")
        return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))

    if request.method == "POST":
        if form.is_valid():
            try:

                handle_uploaded_file(request.FILES['file'])
                input_fasta_path = os.path.join(BASE_DIR, "files", f"{form.cleaned_data['file']}")
                protein_fasta_path = os.path.join(BASE_DIR, "files", f"{request.user.username}.fasta")
                records = SeqIO.parse(input_fasta_path, format="fasta")
                table = form.cleaned_data['translate_table']
                to_stop = form.cleaned_data['to_stop']

                if to_stop is True:

                    protein_file = open(protein_fasta_path, "w")

                    for i in records:
                        seq_dna = i.seq
                        dna_id = i.id
                        dna_description = i.description
                        dna_name = i.name

                        SeqIO.write(
                            SeqRecord(
                                seq=Seq(seq_dna).translate(table=table),
                                id=dna_id.encode().decode(encoding="utf-8", errors="ignore"),
                                description=dna_description,
                                name=dna_name,
                            ),

                            protein_file,

                            'fasta'

                        )
                else:

                    protein_file = open(protein_fasta_path, "w")

                    for i in records:
                        seq_dna = i.seq
                        dna_id = i.id
                        dna_description = i.description
                        dna_name = i.name

                        SeqIO.write(
                            SeqRecord(
                                seq=Seq(seq_dna).translate(table=table).replace("*", ""),
                                id=dna_id.encode().decode(encoding="utf-8", errors="ignore"),
                                description=dna_description,
                                name=dna_name,
                            ),
                            protein_file,

                            'fasta'

                        )

                return HttpResponseRedirect(reverse('bioinformatic:fasta_download'))

            except Bio.BiopythonWarning:
                pass

            finally:
                os.remove(input_fasta_path)

        return redirect('bioinformatic:fasta_protein_download')

    return render(request, "bioinformatic/fasta/read.html", {'form': form, 'bre': "Fasta Dosyası Translasyon"})


files_names = []


def fasta_file_combine(request):
    global writing_fasta
    form = MultipleUploadFileForm(request.POST or None, request.FILES or None)
    if request.method == "POST":
        if form.is_valid():
            files = request.FILES.getlist('file_field')
            combine_fasta = os.path.join(BASE_DIR, "files", "combined.fasta")
            combined_fasta_path = Path(combine_fasta)

            if combined_fasta_path.exists():
                os.remove(combine_fasta)

            for fasta in files:
                handle_uploaded_file(fasta)
                files_names.append(fasta.name)

            path = os.path.join(BASE_DIR, "files", f"{files_names[0]}")

            with open(path, "a") as combine_file:
                for i in files_names[1:]:
                    with open(os.path.join(BASE_DIR, "files", f"{i}"), 'r') as read_fasta:
                        combine_file.write(f"{read_fasta.read()}".replace("\n", ""))
                        read_fasta.close()
                        os.remove(os.path.join(BASE_DIR, "files", f"{i}"))

            os.rename(path, os.path.join(BASE_DIR, "files", "combined.fasta"))

            return redirect('bioinformatic:combine_fasta_download')

        else:
            form = MultipleUploadFileForm()

    return render(request, "bioinformatic/fasta/multifile.html", {'form': form})
