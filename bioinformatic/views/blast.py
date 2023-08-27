import os
from pathlib import Path
from django.contrib.auth.decorators import login_required
from Bio import SearchIO, SeqIO
from django.shortcuts import render, redirect
from bioinformatic.forms.blast import BlastXMLForm

from bioinformatic.models import BlastQueryResults
from Bio.Blast import NCBIWWW, NCBIXML
from django.core.files import File

BASE_DIR = Path(__file__).resolve().parent.parent
path = os.path.join(BASE_DIR, 'files\\')


def handle_uploaded_file(f):
    with open(os.path.join(BASE_DIR, "files", f"{f.name}"), 'wb+') as destination:
        for chunk in f.chunks():
            destination.write(chunk)


@login_required
def fasta_blast_tools(request):
    form = BlastXMLForm(request.POST or None, request.FILES or None)
    if request.method == "POST":
        if form.is_valid():
            handle_uploaded_file(request.FILES['input_file'])
            input_file_path = os.path.join(BASE_DIR, 'files', '{}'.format(form.cleaned_data['input_file']))
            program = form.cleaned_data['program']
            database = form.cleaned_data['database']
            blast_xml_file_path = os.path.join(BASE_DIR, 'files', '{}.xml'.format(program))

            record = next(SeqIO.parse(input_file_path, format="fasta"))

            result_handle = NCBIWWW.qblast(f"{program}", f"{database}", record.format("fasta"))

            save_file = open(blast_xml_file_path, "w")

            save_file.write(result_handle.read())

            save_file.close()

            blast_qresult = SearchIO.read(blast_xml_file_path, "blast-xml")

            blast_records = NCBIXML.parse(open(blast_xml_file_path))

            hsp_file_path = os.path.join(BASE_DIR, 'files', 'hsp.txt')

            file_obj = open(hsp_file_path, "w")

            for hsp in blast_qresult:
                for blast_record in blast_records:
                    for alignment in blast_record.alignments:
                        for hsps in alignment.hsps:
                            file_obj.writelines(f"{hsp}" + 3 * "\n")
                            file_obj.writelines(f"{hsps}" + 3 * "\n")

            file_obj.close()

            if BlastQueryResults.objects.all().filter(user=request.user.id).exists():
                BlastQueryResults.objects.all().filter(user=request.user.id).delete()

            blast_obj = BlastQueryResults()
            blastxml_file_pathlib = Path(blast_xml_file_path)
            blast_hsp_pathlib = Path(hsp_file_path)

            with blastxml_file_pathlib.open('r') as blastxml:
                blast_obj.output_file = File(blastxml, name=blastxml_file_pathlib.name)
                blast_obj.user = request.user
                blast_obj.save()

            with blast_hsp_pathlib.open('r') as hsp_file_obj:
                blast_obj.hsp_file = File(hsp_file_obj, name=blast_hsp_pathlib.name)
                blast_obj.program = program
                blast_obj.save()

            os.remove(input_file_path)
            os.remove(blast_xml_file_path)
            os.remove(hsp_file_path)

            results = BlastQueryResults.objects.all().filter(user=request.user.id).latest('created')
            return render(request, 'bioinformatic/xml/result.html', {'bre': 'Blast Sonuçları', 'results': results})

        else:
            form = BlastXMLForm()

    return render(request, 'bioinformatic/xml/read.html', {'form': form, 'bre': 'Blast Araçları'})
