import os
from django.views import generic
import sys
from Bio.Application import ApplicationError
from pathlib import Path
from django.shortcuts import render, reverse, redirect
from bioinformatic.models import MultipleSequenceAlignment
from Bio import SeqIO, Phylo
from Bio.Align.Applications import ClustalwCommandline
from Bio.Phylo import PhyloXMLIO
from django.core.files import File
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceCalculator
from bioinformatic.forms.alignments import MaximumLikeHoodForm
from Bio.Phylo.PAML import codeml, baseml
from Bio.Phylo.PAML._paml import PamlError
from django.contrib.auth.decorators import login_required
from django.contrib import messages
from django.http import HttpResponseRedirect
import pandas as pd

BASE_DIR = Path(__file__).resolve().parent.parent.parent
path = os.path.join(BASE_DIR, 'bioinformatic', 'files/')


def handle_uploaded_file(f):
    with open(path + f.name, 'wb+') as destination:
        for chunk in f.chunks():
            destination.write(chunk)


class MaximumLikelihood(generic.DetailView):
    template_name = "bioinformatic/alignments/palm_results.html"
    model = MultipleSequenceAlignment


@login_required
def maxlikehood(request, user, method):
    global clustalw2_exe, cml_exe, bml_exe

    if request.user.is_anonymous:
        from django.conf import settings
        messages.error(request, "Lütfen Giriş Yapınız")
        return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))

    form = MaximumLikeHoodForm(request.POST or None, request.FILES or None)

    try:
        obj = MultipleSequenceAlignment.objects.filter(user=request.user).latest('created')

    except MultipleSequenceAlignment.DoesNotExist:
        return redirect('bioinformatic:multiple_sequence_alignments')

    if request.method == "POST":
        if form.is_valid():

            if obj.palm_tools:
                obj.delete()
                return redirect('bioinformatic:multiple_sequence_alignments')

            input_file = os.path.join(BASE_DIR, 'media', 'MultipleSequenceAlignment', f"{request.user}",
                                      f'{request.user}_in_file.fasta')
            output_file = os.path.join(BASE_DIR, 'bioinformatic', 'files', 'aligment.fasta')

            align_file = os.path.join(BASE_DIR, 'bioinformatic', 'files', f'aligned.txt')

            dnd_file = os.path.join(BASE_DIR, 'bioinformatic', 'files', "tree.dnd")
            tree_file = os.path.join(BASE_DIR, 'bioinformatic', 'files', 'tree.xml')
            stats = os.path.join(BASE_DIR, 'bioinformatic', 'files', 'stats.txt')

            scores_path = os.path.join(BASE_DIR, "bioinformatic", "files", "scores.txt")

            max_likelihood_results = os.path.join(BASE_DIR, "bioinformatic", "files", "result_codeml.txt")
            xml_tree_path = os.path.join(BASE_DIR, "media", "MultipleSequenceAlignment", f"{request.user}", 'tree.xml')
            cluster_file_path = os.path.join(BASE_DIR, "media", "MultipleSequenceAlignment", f"{request.user}",
                                             'cluster.csv')

            if sys.platform.startswith('win32'):
                clustalw2_exe = os.path.join(BASE_DIR, 'bioinformatic', 'apps', 'clustalw2.exe')
                cml_exe = os.path.join(BASE_DIR, "bioinformatic", "apps", "paml4.9j", "bin", "codeml.exe")
                bml_exe = os.path.join(BASE_DIR, "bioinformatic", "apps", "paml4.9j", "bin", "baseml.exe")
            elif sys.platform.startswith('linux'):
                clustalw2_exe = os.path.join(BASE_DIR, 'bioinformatic', 'apps', 'clustalw2')
                cml_exe = os.path.join(BASE_DIR, "bioinformatic", "apps", "paml", "bin", "codeml")
                bml_exe = os.path.join(BASE_DIR, "bioinformatic", "apps", "paml", "bin", "baseml")

            if obj.in_file:
                obj.in_file.delete()
            else:
                obj.in_file = File(request.FILES['file'], name="in_file.fasta")
                obj.save()

            molecule_type = form.cleaned_data['molecule_type']
            tree_type = form.cleaned_data['tree_type']
            palm_tools = form.cleaned_data['palm_tools']

            try:

                records = SeqIO.parse(os.path.join(BASE_DIR, "media", "MultipleSequenceAlignment", f"{request.user}",
                                                   f"{request.user}_in_file.fasta"),
                                      "fasta")

                seq_id = []

                for record in records:
                    seq_id.append(record.id)

                if len(seq_id) < 3:
                    return render(request, "bioinformatic/fasta/notfound.html",
                                  {'msg': "Ağaç oluşturmak için en az 3 canlı türü olmalıdır.",
                                   'url': reverse('bioinformatic:multiple_sequence_alignments')})

                clustalw_cline = ClustalwCommandline(
                    clustalw2_exe,
                    infile=input_file,
                    outfile=output_file,
                    align=True,
                    outorder="ALIGNED",
                    convert=True,
                    output="FASTA",
                    newtree=dnd_file,
                    stats=stats
                )

                assert os.path.isfile(clustalw2_exe), "Clustal W executable missing"
                stdout, stderr = clustalw_cline()

                AlignIO.convert(output_file, 'fasta', align_file, "clustal", molecule_type=molecule_type)
                alignment = AlignIO.read(align_file, 'clustal')

                calculator = DistanceCalculator('identity')
                constructor = DistanceTreeConstructor(calculator, method=tree_type)
                tree = constructor.build_tree(alignment)

                Phylo.write(tree, xml_tree_path, "phyloxml")

                read_stats = open(stats, 'r').readlines()[1:16]

                os.remove(stats)

                for i in read_stats:
                    open(stats, 'a').writelines(i)

                open(scores_path, "w").writelines(stdout)

                scores = open(scores_path, 'r').readlines()[44:52]

                os.remove(scores_path)

                for i in scores:
                    open(scores_path, 'a').writelines(i)

                xml_file_read = PhyloXMLIO.parse(xml_tree_path)

                cols = ["organizma", "uzunluk"]
                rows = []

                for i in xml_file_read:
                    for j in i.clade.find_clades():
                        rows.append({
                            'organizma': j.name,
                            'uzunluk': j.branch_length
                        })

                df = pd.DataFrame(rows, columns=cols)
                df.to_csv(cluster_file_path)

                rows.clear()

                obj.out_file = File(Path(output_file).open('r'), name="out_alignment.fasta")
                obj.tree_type = tree_type
                obj.cluster_csv = File(Path(cluster_file_path).open('r'), name="cluster.csv")
                obj.align_file = File(Path(align_file).open('r'), name="aligned.txt")
                obj.molecule_type = molecule_type
                obj.palm_tools = palm_tools
                obj.stats = File(Path(stats).open('r'), name="stats.txt")
                obj.scores = File(Path(scores_path).open('r'), name="scores.txt")
                obj.tree_file = File(Path(xml_tree_path).open('r'), name="tree.xml")
                obj.save()

                if palm_tools == "codeml":
                    try:
                        cml = codeml.Codeml()
                        cml.alignment = os.path.join(BASE_DIR, "bioinformatic", "files", "aligment.fasta")
                        cml.tree = os.path.join(BASE_DIR, "bioinformatic", "files", "tree.dnd")
                        cml.ctl_file = os.path.join(BASE_DIR, "bioinformatic", "files", f"{palm_tools}.ctl")
                        cml.out_file = os.path.join(BASE_DIR, "bioinformatic", "files", "result_codeml.txt")
                        cml.working_dir = os.path.join(BASE_DIR, "bioinformatic", "files")
                        cml.run(command=cml_exe, verbose=True)

                    except PamlError:
                        pass

                    finally:
                        os.remove(os.path.join(BASE_DIR, "bioinformatic", "files", "lnf"))
                        os.remove(os.path.join(BASE_DIR, "bioinformatic", "files", "rst"))
                        os.remove(os.path.join(BASE_DIR, "bioinformatic", "files", "rub"))
                        os.remove(os.path.join(BASE_DIR, "bioinformatic", "files", "rst1"))

                elif palm_tools == "baseml":
                    try:
                        bml = baseml.Baseml()
                        bml.alignment = os.path.join(BASE_DIR, "bioinformatic", "files", "aligment.fasta")
                        bml.tree = os.path.join(BASE_DIR, "bioinformatic", "files", "tree.dnd")
                        bml.ctl_file = os.path.join(BASE_DIR, "bioinformatic", "files", f"{palm_tools}.ctl")
                        bml.out_file = os.path.join(BASE_DIR, "bioinformatic", "files", "result_codeml.txt")
                        bml.working_dir = os.path.join(BASE_DIR, "bioinformatic", "files")
                        bml.run(command=bml_exe, verbose=True)

                    except PamlError:
                        pass

                    finally:
                        os.remove(os.path.join(BASE_DIR, "bioinformatic", "files", "2base.t"))
                        os.remove(os.path.join(BASE_DIR, "bioinformatic", "files", "lnf"))

                with Path(max_likelihood_results).open(mode='r') as max_file:
                    obj.ml_file = File(max_file, name=f"results_{palm_tools}.txt")
                    obj.save()

                results = MultipleSequenceAlignment.objects.all().filter(user=request.user).latest('created')

                os.remove(input_file)
                os.remove(output_file)
                os.remove(stats)
                os.remove(scores_path)
                os.remove(align_file)
                os.remove(dnd_file)
                os.remove(xml_tree_path)
                os.remove(max_likelihood_results)
                os.remove(os.path.join(BASE_DIR, "media", "MultipleSequenceAlignment", "codeml.ctl"))
                os.remove(os.path.join(BASE_DIR, "bioinformatic", "files", "rst"))
                os.remove(os.path.join(BASE_DIR, "bioinformatic", "files", "rub"))
                os.remove(os.path.join(BASE_DIR, "bioinformatic", "files", "rst1"))

                return HttpResponseRedirect(
                    reverse('bioinformatic:msa_results',
                            args=(request.user, obj.method, obj.molecule_type.lower(), obj.pk)))

            except ApplicationError:
                os.remove(
                    os.path.join(BASE_DIR, 'bioinformatic', 'files', '{}'.format(form.cleaned_data['file'])))
                os.remove(os.path.join(BASE_DIR, 'bioinformatic', 'files', 'aligned.fasta'))

                return render(request, 'bioinformatic/fasta/notfound.html', {
                    'msg': 'Hatalı Dosya Seçtiniz. Lütfen fasta dosyası seçiniz.',
                    'url': reverse('bioinformatic:multiple_sequence_alignments')})

    return render(request, "bioinformatic/alignments/palm.html", {'form': form, 'bre': 'Maximum Likelihood '})
