from django.shortcuts import render, redirect, reverse
from bioinformatic.forms.file import FileReadForm
from bioinformatic.forms.writing import GenbankWritingForm
from bioinformatic.forms.add import AddFastaData
from Bio import SeqIO
from bioinformatic.models import GenbankRead
from pathlib import Path
import os
from django.views import generic
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import plotly.express as px

from bioinformatic.forms.translation import GenbankTranslationForm

BASE_DIR = Path(__file__).resolve().parent.parent
path = os.path.join(BASE_DIR, 'files\\')


def handle_uploaded_file(f):
    with open(path + f.name, 'wb+') as destination:
        for chunk in f.chunks():
            destination.write(chunk)


def genbank_read(request):
    global file, feature, records
    form1 = FileReadForm(request.POST or None, request.FILES or None)
    form2 = GenbankTranslationForm(request.POST or None)
    if request.method == "POST":

        if form1.is_valid():

            handle_uploaded_file(request.FILES['file'])

            file = os.path.join(BASE_DIR, 'files\\{}'.format(form1.cleaned_data['file']))

            try:

                file_obj = open(file, 'r')

                records = SeqIO.parse(file, "genbank")

                file_obj.close()

                if GenbankRead.objects.exists():
                    GenbankRead.objects.all().delete()

                    for record in records:
                        for feature in record.features:
                            if feature.type == "CDS":
                                if feature.qualifiers.get('translation') is not None:
                                    GenbankRead.objects.create(
                                        protein_sequence=feature.qualifiers.get('translation')[0],
                                        protein_id=feature.qualifiers.get('protein_id')[0],
                                        protein_sequence_len=len(feature.qualifiers.get('translation')[0]),
                                        dna_sequence=record.seq,
                                        dna_sequence_len=len(record.seq),
                                        organism=record.annotations['organism'],
                                        taxonomy=record.annotations['taxonomy'],
                                        description=record.description,
                                    )

                            else:
                                GenbankRead.objects.create(
                                    organism=record.annotations['organism'],
                                    taxonomy=record.annotations['taxonomy'],
                                    description=record.description,
                                    dna_sequence=record.seq,
                                    dna_sequence_len=len(record.seq),
                                )

                else:
                    for record in records:
                        for feature in record.features:
                            if feature.type == "CDS":
                                if feature.qualifiers.get('translation') is not None:
                                    GenbankRead.objects.create(
                                        protein_sequence=feature.qualifiers.get('translation')[0],
                                        protein_id=feature.qualifiers.get('protein_id')[0],
                                        protein_sequence_len=len(feature.qualifiers.get('translation')[0]),
                                        dna_sequence=record.seq,
                                        dna_sequence_len=len(record.seq),
                                        organism=record.annotations['organism'],
                                        taxonomy=record.annotations['taxonomy'],
                                        description=record.description,
                                    )

                            else:
                                GenbankRead.objects.create(
                                    organism=record.annotations['organism'],
                                    taxonomy=record.annotations['taxonomy'],
                                    description=record.description,
                                    dna_sequence=record.seq,
                                    dna_sequence_len=len(record.seq),
                                )

            except UnicodeDecodeError:
                os.remove(file)
                return render(request, 'bioinformatic/fasta/notfound.html',
                              {'msg': 'Hatalı Dosya Seçtiniz!',
                               'url': reverse('bioinformatic:genbank_read')})

            finally:
                open(file, 'r').close()
                os.remove(file)

            if GenbankRead.objects.exists():
                return redirect('bioinformatic:genbank_region')
            else:
                return render(request, 'bioinformatic/genbank/notfound.html',
                              {'msg': 'Hatalı Dosya. Lütfen Genbank Dosyası Giriniz.',
                               'url': reverse('bioinformatic:genbank_read')})

        else:

            form1 = FileReadForm(request.POST or None, request.FILES or None)
            form2 = GenbankTranslationForm(request.POST or None)

            os.remove(file)

            return redirect('bioinformatic:genbank_read')

    return render(request, 'bioinformatic/genbank/read.html',
                  {'form1': form1, 'bre': 'Genbank Dosyası Okuması'})


def delete_genbank(request):
    GenbankRead.objects.all().delete()
    return redirect('bioinformatic:genbank_read')


class GenBankResultView(generic.ListView):
    template_name = "bioinformatic/genbank/list.html"
    model = GenbankRead
    paginate_by = 100

    def get_context_data(self, *, object_list=None, **kwargs):
        context = super(GenBankResultView, self).get_context_data(**kwargs)

        organism = []

        dna_sequence_len = []

        protein_id = []

        protein_seq_len = []

        data = {
            "organism": organism,
            "dna_sequence_len": dna_sequence_len,
            "protein_id": protein_id
        }

        for seq in GenbankRead.objects.all():
            dna_sequence_len.append(seq.dna_sequence_len)

        for organizma in GenbankRead.objects.all():
            organism.append(organizma.organism)

        for protein in GenbankRead.objects.all():
            protein_id.append(protein.protein_id)

        for protein_seq in GenbankRead.objects.all():
            protein_seq_len.append(protein_seq.protein_sequence_len)

        fig = px.bar(x=organism, y=dna_sequence_len, category_orders={'organism': organism}, color=organism,
                     title="GenBank Dosyası")

        context['count'] = GenbankRead.objects.all().count()
        context['bre'] = "Genbank Dosya Okuması Sonuçları"
        context['fig'] = fig

        return context


class GenbankDetailView(generic.DetailView):
    model = GenbankRead
    template_name = "bioinformatic/genbank/detail.html"

    def get_context_data(self, **kwargs):
        context = super(GenbankDetailView, self).get_context_data(**kwargs)
        context['bre'] = 'Genbank okuması'
        return context


def genbank_writing(request):
    form = GenbankWritingForm(request.POST or None)
    if request.method == "POST":
        if form.is_valid():

            lokus = form.cleaned_data["lokus"]
            id = form.cleaned_data["id"]
            descriptions = form.cleaned_data["description"]
            sequence = form.cleaned_data["sequence"]
            molecule_type = form.cleaned_data["molecule_type"]
            keywords = form.cleaned_data["keywords"]
            taxonomy = form.cleaned_data["taxonomy"]
            references = form.cleaned_data["references"]
            dbxref = form.cleaned_data["dbxref"]
            source = form.cleaned_data["source"]
            organism = form.cleaned_data["organism"]
            features = form.cleaned_data["features"]
            sequence = Seq(sequence).upper()

            bad_chars = [';', ':', '!', "*", "\n", '"', "\r", " "]

            for i in bad_chars:
                sequence = sequence.replace(i, '')

            rec2 = SeqRecord(
                sequence,
                id=id,
                description=descriptions, name=lokus, dbxrefs=[].append(dbxref), features=[].append(features)
            )

            rec2.annotations["molecule_type"] = molecule_type
            rec2.annotations["keywords"] = keywords
            rec2.annotations["taxonomy"] = taxonomy
            rec2.annotations["lokus"] = lokus
            rec2.annotations["source"] = source
            rec2.annotations["organism"] = organism
            rec2.annotations["references"] = references

            file = os.path.join(BASE_DIR, 'files\\file.gbk')

            SeqIO.write(rec2, file, "genbank")

            return redirect("bioinformatic:genbank_download")

        else:
            msg = "Bir hata meydana geldi"

            return render(request,
                          "bioinformatic/fasta/notfound.html", {
                              "msg": msg
                          })

    return render(request, "bioinformatic/genbank/writing.html", {
        "form": form,
        "bre": "Genbank Dosyası Yazma"
    })


def append_multiple_lines(file_name, lines_to_append):
    # Open the file in append & read mode ('a+')
    with open(file_name, "a+") as file_object:
        appendEOL = False
        # Move read cursor to the start of file.
        file_object.seek(0)
        # Check if file is not empty
        data = file_object.read(100)
        if len(data) > 0:
            appendEOL = True
        # Iterate over each string in the list
        for line in lines_to_append:
            # If file is not empty then append '\n' before first line for
            # other lines always append '\n' before appending line
            if appendEOL == True:
                file_object.write("\n")
            else:
                appendEOL = True
            # Append element at the end of file
            file_object.write(line)
            file_object.close()


def fasta_add(request):
    form = AddFastaData(request.POST, request.FILES)
    if request.method == "POST":
        if form.is_valid():
            try:

                handle_uploaded_file(request.FILES["file"])

                id = form.cleaned_data["id"]
                descriptions = form.cleaned_data["description"]
                sequence = form.cleaned_data["sequence"]
                sequence = Seq(sequence)

                bad_chars = [';', ':', '!', "*", "\n", '"', "\r"]

                for i in bad_chars:
                    sequence = sequence.replace(i, '')

                rec1 = SeqRecord(sequence, id=id, description=descriptions)

                print(rec1)

                file = os.path.join(BASE_DIR, 'files\\{}'.format(request.FILES['file']))
                new = os.path.join(BASE_DIR, 'files\\file.fasta')

                os.rename(file, new)

                append_multiple_lines(new, rec1)

                return redirect("bioinformatic:download")

            except:
                pass
        else:

            msg = "Bir hata meydana geldi"

            return render(request, "bioinformatic/fasta/notfound.html", {
                "msg": msg
            })

    return render(request, "bioinformatic/fasta/add.html", {
        "form": form,
        "bre": "Fasta Dosyası Veri Ekleme"
    })
