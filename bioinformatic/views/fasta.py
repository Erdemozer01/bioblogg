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
from bioinformatic.models import GraphicModels
from django.core.files import File
from django.contrib.auth.decorators import login_required
from django.views import generic
from django.contrib import messages

BASE_DIR = Path(__file__).resolve().parent.parent
path = os.path.join(BASE_DIR, 'files\\')


def handle_uploaded_file(f):
    with open(os.path.join(BASE_DIR, "files", f"{f}"), 'wb+') as destination:
        for chunk in f.chunks():
            destination.write(chunk)


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
    if request.user.is_anonymous:
        from django.conf import settings
        messages.error(request, "Lütfen Giriş Yapınız")
        return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))
    fastaform = FastaWritingForm(request.POST or None)
    if request.method == "POST":
        if fastaform.is_valid():

            id = fastaform.cleaned_data["id"]
            descriptions = fastaform.cleaned_data["description"]
            sequence = fastaform.cleaned_data["sequence"]
            sequence = Seq(sequence)
            bad_chars = [';', ':', '!', "*", "\n", '"', "\r"]

            for i in bad_chars:
                sequence = sequence.replace(i, '')

            rec1 = SeqRecord(
                sequence,
                id=id,
                description=descriptions
            )

            file = os.path.join(BASE_DIR, f'files\\{request.user.username}.fasta')

            SeqIO.write(rec1, file, "fasta")

            return redirect("bioinformatic:fasta_download")

        else:

            msg = "Bir hata meydana geldi"

            return render(request, 'bioinformatic/fasta/notfound.html', {
                "msg": msg
            })

    return render(request, "bioinformatic/fasta/writing.html", {
        "form": fastaform,
        "bre": "Fasta Dosyası Yazma"
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
    form = DNAFastaFileTranslateForm(request.POST or None, request.FILES or None)
    if request.method == "POST":
        if form.is_valid():
            try:
                handle_uploaded_file(request.FILES['file'])
                input_fasta_path = os.path.join(BASE_DIR, "files", f"{form.cleaned_data['file']}")
                protein_fasta_path = os.path.join(BASE_DIR, "files", "protein.txt")
                records = SeqIO.parse(input_fasta_path, format="fasta")
                table = form.cleaned_data['translate_table']

                protein_file = open(protein_fasta_path, "w")

                for record in records:
                    protein_file.write(f"{record.description}")
                    protein_file.write("\n")
                    protein_file.write("Sekans Uzunlugu: " + f"{len(record.translate().seq)}")
                    protein_file.write("\n")
                    protein_file.write(f"{record.translate(table=table).seq}")
                    protein_file.write(2 * "\n")

                os.remove(input_fasta_path)

            except Bio.BiopythonWarning:
                pass

        return redirect('bioinformatic:fasta_protein_download')

    return render(request, "bioinformatic/fasta/translate.html", {'form': form, 'bre':"Fasta Dosyası Translasyon"})


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
