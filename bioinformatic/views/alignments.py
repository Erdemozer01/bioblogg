import os.path
import sys
from django.contrib import messages
from django.http import HttpResponseRedirect
from django.core.files import File
from pathlib import Path
import Bio.Application
from Bio import pairwise2, Phylo, SeqIO
from Bio.Align import substitution_matrices
from django.shortcuts import render, redirect, reverse
from Bio import AlignIO
import subprocess
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceCalculator
from Bio.Align.Applications import ClustalwCommandline, ClustalOmegaCommandline, MuscleCommandline
from bioinformatic \
    .forms.alignments import GlobalForm, LocalForm, MultipleSequenceAlignmentForm, MultipleSequenceAlignmentSelectForm, \
    MaximumLikeHoodForm
from bioinformatic.models import MultipleSequenceAlignment, MaximumFileSize
from django.views import generic
from dash import html, dcc
import dash_bio as dashbio
from dash.dependencies import Input, Output
from django_plotly_dash import DjangoDash
import json, math
import dash_cytoscape as cyto
import pandas as pd
from Bio.Phylo import PhyloXMLIO
from django.urls import reverse_lazy

BASE_DIR = Path(__file__).resolve().parent.parent.parent


def handle_uploaded_file(f):
    with open(os.path.join(BASE_DIR, "bioinformatic", "files", f"{f}"), 'wb+') as destination:
        for chunk in f.chunks():
            destination.write(chunk)


def global_alignment(request):
    if MaximumFileSize.objects.exists():
        max_file_size = MaximumFileSize.objects.all().latest('created')
    form = GlobalForm(request.POST or None)

    if request.POST:

        if form.is_valid():

            seq1 = form.cleaned_data["seq1"]
            seq2 = form.cleaned_data["seq2"]
            bad_chars = [';', ':', '!', "*", "\n", '"', "\r", " "]
            for i in bad_chars:
                seq1 = seq1.replace(i, '')
                seq2 = seq2.replace(i, '')

            alignments = pairwise2.align.globalxx(seq1, seq2)
            result = pairwise2.format_alignment(*alignments[0], full_sequences=True)
            from Bio import Align
            aligner = Align.PairwiseAligner()
            alignments2 = aligner.align(seq1, seq2)

            filepath = os.path.join(BASE_DIR, 'bioinformatic\\files\\global_alignment.txt')
            file = open(filepath, "w")

            file.write(result)
            file.write("\n\n")
            file.write(str(alignments2[0].substitutions))

            file.close()

            return redirect("bioinformatic:global_alignment_download")

        else:
            form = GlobalForm()

    return render(request, "bioinformatic/alignments/global.html", {"form": form, "bre": "Global Alignment"})


def local_alignment(request):
    form = LocalForm(request.POST or None)
    if request.POST:
        if form.is_valid():
            seq1 = form.cleaned_data["seq1"]
            seq2 = form.cleaned_data["seq2"]
            align = form.cleaned_data["algo"]

            bad_chars = [';', ':', '!', "*", "\n", '"', "\r", " "]
            for i in bad_chars:
                seq1 = seq1.replace(i, '')
                seq2 = seq2.replace(i, '')

            matrice = substitution_matrices.load(f"{align}")
            alignments = pairwise2.align.localds(seq1, seq2, matrice, -10, -1)
            result = pairwise2.format_alignment(*alignments[0], full_sequences=True)

            filepath = os.path.join(BASE_DIR, 'bioinformatic\\files\\local_alignment.txt')
            file = open(filepath, "w")

            file.write(f"{align} Matrisine göre dizileme yapılmıştır\n\n")
            file.close()

            app = open(filepath, "a")
            app.write(result)

            app.close()
            return redirect("bioinformatic:local_alignment_download")

        else:
            form = LocalForm()

    return render(request, "bioinformatic/alignments/local.html", {"form": form, "bre": "Local Alignment"})


def MultipleSequenceAlignmentView(request):
    if request.user.is_anonymous:
        from django.conf import settings
        messages.error(request, "Lütfen giriş yapınız")
        return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))

    obj = MultipleSequenceAlignment()
    form = MultipleSequenceAlignmentSelectForm(request.POST or None)
    if request.method == "POST":
        if form.is_valid():
            if MultipleSequenceAlignment.objects.all().filter(user=request.user).exists():
                MultipleSequenceAlignment.objects.all().filter(user=request.user).all().delete()
            method = form.cleaned_data['method']
            obj.user = request.user
            obj.method = method
            obj.save()
            if obj.method == "paml":
                return HttpResponseRedirect(
                    reverse('bioinformatic:maximum_likelihood',
                            args=(request.user, obj.method)))
            else:
                return HttpResponseRedirect(
                    reverse('bioinformatic:multiple_sequence_alignments_analiz',
                            args=(request.user, obj.method)))

    return render(request, "bioinformatic/alignments/msa_select.html", {'form': form})


def MultipleSeqAlignment(request, user, method):
    global clustalw2_exe, muscle_exe, clustal_omega_exe, clustal_result, clustalw2, clustalx_exe, cml_exe,obj

    if request.user.is_anonymous:
        from django.conf import settings
        messages.error(request, "Lütfen Giriş Yapınız")
        return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))

    form = MultipleSequenceAlignmentForm(request.POST or None, request.FILES or None)

    try:
        obj = MultipleSequenceAlignment.objects.filter(user=request.user).latest('created')

    except MultipleSequenceAlignment.DoesNotExist:
        return redirect('bioinformatic:multiple_sequence_alignments')

    if request.method == "POST":
        if form.is_valid():
            if obj:
               obj.delete()

            obj.in_file = File(request.FILES['file'], name="in_file.fasta")
            obj.save()

            molecule_type = form.cleaned_data['molecule_type']
            tree_type = form.cleaned_data['tree_type']
            alignment_filetype = form.cleaned_data['alignment_filetype']

            in_file_path = os.path.join(BASE_DIR, "media", "MultipleSequenceAlignment",
                                        f"{request.user}", f"{request.user}_in_file.fasta")

            out_file_path = os.path.join(BASE_DIR, "media", "MultipleSequenceAlignment", f"{request.user}",
                                         'alignment.fasta')
            align_file_path = os.path.join(BASE_DIR, "media", "MultipleSequenceAlignment", f"{request.user}",
                                           f'aligned.{alignment_filetype}')
            stats_file_path = os.path.join(BASE_DIR, "media", "MultipleSequenceAlignment", f"{request.user}",
                                           'stats.txt')
            scores_path = os.path.join(BASE_DIR, "media", "MultipleSequenceAlignment", f"{request.user}",
                                       "scores.txt")
            xml_tree_path = os.path.join(BASE_DIR, "media", "MultipleSequenceAlignment", f"{request.user}",
                                         'tree.xml')
            cluster_file_path = os.path.join(BASE_DIR, "media", "MultipleSequenceAlignment", f"{request.user}",
                                             'cluster.csv')

            if alignment_filetype == "phylip":
                align_file_path = os.path.join(BASE_DIR, "media", "MultipleSequenceAlignment",
                                               f"{request.user}",
                                               'aligned.phy')

            if obj.method == "muscle":

                try:
                    if sys.platform.startswith('win32'):
                        muscle_exe = os.path.join(BASE_DIR, 'bioinformatic', 'apps', 'muscle3.8.425_win32.exe')
                        muscle_cline = MuscleCommandline(
                            muscle_exe,
                            input=in_file_path,
                            out=out_file_path,
                        )
                        stdin, stdout = muscle_cline()
                    elif sys.platform.startswith('linux'):
                        muscle_exe = os.path.join(BASE_DIR, 'bioinformatic', 'apps', 'muscle3.8.425_i86linux32')
                        muscle_result = subprocess.check_output(
                            [muscle_exe, "-in", in_file_path, "-out", out_file_path])

                    AlignIO.convert(out_file_path, 'fasta', align_file_path, f'{alignment_filetype}',
                                    molecule_type=molecule_type)

                    alignment = AlignIO.read(align_file_path, f'{alignment_filetype}')

                    calculator = DistanceCalculator('identity')

                    constructor = DistanceTreeConstructor(calculator, method=tree_type)
                    trees = constructor.build_tree(alignment)

                    Phylo.write(trees, xml_tree_path, "phyloxml")

                    records = SeqIO.parse(in_file_path, "fasta")

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

                    seq_id = []

                    for record in records:
                        seq_id.append(record.id)

                    if len(seq_id) < 3:
                        return render(request, "bioinformatic/fasta/notfound.html",
                                      {'msg': "Ağaç oluşturmak için en az 3 canlı türü olmalıdır.",
                                       'url': reverse('bioinformatic:multiple_sequence_alignments')})

                    obj.tree_type = tree_type
                    obj.alignment_filetype = alignment_filetype
                    obj.molecule_type = molecule_type
                    obj.align_file = File(Path(align_file_path).open('r'),
                                          name="aligned.{}".format(alignment_filetype))
                    obj.out_file = File(Path(out_file_path).open('r'), name="out_alignment.fasta")
                    obj.tree_file = File(Path(xml_tree_path).open('r'), name="tree.xml")
                    obj.cluster_csv = File(Path(cluster_file_path).open('r'), name="cluster.csv")
                    obj.save()

                    handle = open(in_file_path)
                    handle.close()
                    os.remove(in_file_path)
                    handle = open(out_file_path)
                    handle.close()
                    os.remove(out_file_path)
                    handle = open(xml_tree_path)
                    handle.close()
                    os.remove(xml_tree_path)
                    handle = open(align_file_path)
                    handle.close()
                    os.remove(align_file_path)
                    handle = open(cluster_file_path)
                    handle.close()
                    os.remove(cluster_file_path)

                except Bio.Application.ApplicationError:
                    try:
                        os.remove(
                            os.path.join(BASE_DIR, 'bioinformatic', 'files',
                                         '{}'.format(form.cleaned_data['file'])))
                        os.remove(os.path.join(BASE_DIR, 'bioinformatic', 'files', 'aligned.fasta'))

                        return render(request, 'bioinformatic/fasta/notfound.html', {
                            'msg': 'Hatalı Dosya Seçtiniz. Lütfen fasta dosyası seçiniz.',
                            'url': reverse('bioinformatic:multiple_sequence_alignments')})

                    except FileNotFoundError:
                        return render(request, 'bioinformatic/fasta/notfound.html', {
                            'msg': 'Dosya Bulunamadı',
                            'url': reverse('bioinformatic:multiple_sequence_alignments')})

            elif obj.method == "clustalw2":
                try:
                    if sys.platform.startswith('win32'):
                        clustalw2_exe = os.path.join(BASE_DIR, 'bioinformatic', 'apps', 'clustalw2.exe')

                    elif sys.platform.startswith('linux'):
                        clustalw2_exe = os.path.join(BASE_DIR, 'bioinformatic', 'apps', 'clustalw2')

                    clustalw_cline = ClustalwCommandline(
                        clustalw2_exe,
                        infile=in_file_path,
                        outfile=out_file_path,
                        align=True,
                        outorder="ALIGNED",
                        convert=True,
                        output="FASTA",
                        stats=stats_file_path,
                    )

                    assert os.path.isfile(clustalw2_exe), "Clustal W executable missing"
                    stdout, stderr = clustalw_cline()

                    AlignIO.convert(out_file_path, 'fasta', align_file_path, f'{alignment_filetype}',
                                    molecule_type=molecule_type)
                    alignment = AlignIO.read(align_file_path, f'{alignment_filetype}')

                    calculator = DistanceCalculator('identity')

                    constructor = DistanceTreeConstructor(calculator, method=tree_type)
                    tree = constructor.build_tree(alignment)

                    Phylo.write(tree, xml_tree_path, "phyloxml")

                    records = SeqIO.parse(in_file_path, "fasta")

                    seq_id = []

                    for record in records:
                        seq_id.append(record.id)

                    if len(seq_id) < 3:
                        return render(request, "bioinformatic/fasta/notfound.html",
                                      {'msg': "Ağaç oluşturmak için en az 3 canlı türü olmalıdır.",
                                       'url': reverse('bioinformatic:multiple_sequence_alignments')})

                    open(scores_path, "w").writelines(stdout)

                    scores = open(scores_path, 'r').readlines()[44:52]

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

                    obj.user = request.user
                    obj.method = method
                    obj.tree_type = tree_type
                    obj.alignment_filetype = alignment_filetype
                    obj.molecule_type = molecule_type
                    obj.align_file = File(Path(align_file_path).open('r'), name="aligned.phy")
                    obj.out_file = File(Path(out_file_path).open('r'), name="out_alignment.fasta")
                    obj.tree_file = File(Path(xml_tree_path).open('r'), name="tree.xml")
                    obj.cluster_csv = File(Path(cluster_file_path).open('r'), name="cluster.csv")
                    obj.stats = File(Path(stats_file_path).open('r'), name="stats.txt")
                    obj.scores = File(Path(scores_path).open('r'), name="scores.txt")
                    obj.save()

                    handle = open(in_file_path)
                    handle.close()
                    os.remove(in_file_path)
                    handle = open(out_file_path)
                    handle.close()
                    os.remove(out_file_path)
                    handle = open(align_file_path)
                    handle.close()
                    os.remove(align_file_path)
                    handle = open(xml_tree_path)
                    handle.close()
                    os.remove(xml_tree_path)
                    handle = open(stats_file_path)
                    handle.close()
                    os.remove(stats_file_path)
                    handle = open(scores_path)
                    handle.close()
                    os.remove(scores_path)

                except Bio.Application.ApplicationError:
                    obj.delete()

                    return render(request, 'bioinformatic/fasta/notfound.html', {
                        'msg': 'Hatalı Dosya Seçtiniz. Lütfen fasta dosyası seçiniz.',
                        'url': reverse('bioinformatic:multiple_sequence_alignments')})

            elif obj.method == "omega":
                try:
                    if sys.platform.startswith('win32'):
                        clustal_omega_exe = os.path.join(BASE_DIR, 'bioinformatic', 'apps',
                                                         'clustal-omega-1.2.2-win64/clustalo.exe')
                    elif sys.platform.startswith('linux'):
                        clustal_omega_exe = os.path.join(BASE_DIR, 'bioinformatic', 'apps',
                                                         'clustalo-1.2.4-Ubuntu-32-bit')

                    clustal_omega_cline = ClustalOmegaCommandline(
                        clustal_omega_exe,
                        infile=in_file_path,
                        outfile=out_file_path,
                        force=True,
                        verbose=True,
                        auto=True,
                        usekimura="yes"
                    )

                    if sys.platform.startswith('win32'):
                        assert os.path.isfile(clustal_omega_exe)
                        stdout, stderr = clustal_omega_cline()

                    elif sys.platform.startswith('linux'):
                        subprocess.Popen(str(clustal_omega_cline), stdin=subprocess.PIPE,
                                         stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE, universal_newlines=True,
                                         shell=(sys.platform != "win32"))

                    AlignIO.convert(out_file_path, 'fasta', align_file_path, f'{alignment_filetype}',
                                    molecule_type=molecule_type)
                    alignment = AlignIO.read(align_file_path, f'{alignment_filetype}')

                    calculator = DistanceCalculator('identity')

                    constructor = DistanceTreeConstructor(calculator, method=tree_type)

                    tree = constructor.build_tree(alignment)

                    Phylo.write(tree, xml_tree_path, "phyloxml")

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

                    records = SeqIO.parse(in_file_path, "fasta")

                    seq_id = []

                    for record in records:
                        seq_id.append(record.id)

                    if len(seq_id) < 3:
                        return render(request, "bioinformatic/fasta/notfound.html",
                                      {'msg': "Ağaç oluşturmak için en az 3 canlı türü olmalıdır.",
                                       'url': reverse('bioinformatic:multiple_sequence_alignments')})

                    obj.tree_type = tree_type
                    obj.alignment_filetype = alignment_filetype
                    obj.molecule_type = molecule_type
                    obj.align_file = File(Path(align_file_path).open('r'),
                                          name="aligned.{}".format(alignment_filetype))
                    obj.out_file = File(Path(out_file_path).open('r'), name="out_alignment.fasta")
                    obj.tree_file = File(Path(xml_tree_path).open('r'), name="tree.xml")
                    obj.cluster_csv = File(Path(cluster_file_path).open('r'), name="cluster.csv")
                    obj.save()

                    handle = open(in_file_path)
                    handle.close()
                    os.remove(in_file_path)
                    handle = open(out_file_path)
                    handle.close()
                    os.remove(out_file_path)
                    handle = open(align_file_path)
                    handle.close()
                    os.remove(align_file_path)
                    handle = open(xml_tree_path)
                    handle.close()
                    os.remove(xml_tree_path)

                except Bio.Application.ApplicationError:
                    os.remove(
                        os.path.join(BASE_DIR, 'bioinformatic', 'files',
                                     '{}'.format(form.cleaned_data['file'])))
                    os.remove(os.path.join(BASE_DIR, 'bioinformatic', 'files', 'aligned.fasta'))

                    return render(request, 'bioinformatic/fasta/notfound.html', {
                        'msg': 'Hatalı Dosya Seçtiniz. Lütfen fasta dosyası seçiniz.',
                        'url': reverse('bioinformatic:multiple_sequence_alignments')})

            return HttpResponseRedirect(
                reverse('bioinformatic:msa_results',
                        args=(request.user, obj.method, obj.molecule_type.lower(), obj.pk)))

        else:
            form = MultipleSequenceAlignmentForm()

    return render(request, 'bioinformatic/alignments/msa_analiz.html',
                  {'form': form, 'bre': 'Multiple Sekans Alignment', 'obj': obj})


class MultipleSeqDetailView(generic.DetailView):
    template_name = "bioinformatic/alignments/msa_results.html"
    model = MultipleSequenceAlignment

    def get_context_data(self, **kwargs):
        context = super(MultipleSeqDetailView, self).get_context_data(**kwargs)
        context['bre'] = "Multiple Sekans Alignment Sonuçları"
        return context


class PhyloTreeDetailView(generic.DetailView):
    template_name = "bioinformatic/alignments/tree.html"
    model = MultipleSequenceAlignment

    def get_queryset(self):
        return MultipleSequenceAlignment.objects.filter(user=self.request.user)

    def get(self, request, *args, **kwargs):
        global tree
        if self.request.user.is_anonymous:
            from django.contrib import messages
            messages.error(request, "Lütfen Giriş Yapınız")
            return redirect('login')
        else:
            pass

        app = DjangoDash('PhylogeneticTree')
        tree_file = os.path.join(BASE_DIR, "media", "MultipleSequenceAlignment", f"{self.request.user}",
                                 f"{self.request.user}_tree.xml")

        def generate_elements(tree, xlen=30, ylen=30, grabbable=False):
            def get_col_positions(tree, column_width=80):
                taxa = tree.get_terminals()

                # Some constants for the drawing calculations
                max_label_width = max(len(str(taxon)) for taxon in taxa)
                drawing_width = column_width - max_label_width - 1

                """Create a mapping of each clade to its column position."""
                depths = tree.depths()
                # If there are no branch lengths, assume unit branch lengths
                if not max(depths.values()):
                    depths = tree.depths(unit_branch_lengths=True)
                    # Potential drawing overflow due to rounding -- 1 char per tree layer
                fudge_margin = int(math.ceil(math.log(len(taxa), 2)))
                cols_per_branch_unit = ((drawing_width - fudge_margin) /
                                        float(max(depths.values())))
                return dict((clade, int(blen * cols_per_branch_unit + 1.0))
                            for clade, blen in depths.items())

            def get_row_positions(tree):
                taxa = tree.get_terminals()
                positions = dict((taxon, 2 * idx) for idx, taxon in enumerate(taxa))

                def calc_row(clade):
                    for subclade in clade:
                        if subclade not in positions:
                            calc_row(subclade)
                    positions[clade] = ((positions[clade.clades[0]] +
                                         positions[clade.clades[-1]]) // 2)

                calc_row(tree.root)
                return positions

            def add_to_elements(clade, clade_id):
                children = clade.clades

                pos_x = col_positions[clade] * xlen
                pos_y = row_positions[clade] * ylen

                cy_source = {
                    "data": {"id": clade_id},
                    'position': {'x': pos_x, 'y': pos_y},
                    'classes': 'nonterminal',
                    'grabbable': grabbable
                }
                nodes.append(cy_source)

                if clade.is_terminal():
                    cy_source['data']['name'] = clade.name
                    cy_source['classes'] = 'terminal'

                for n, child in enumerate(children):
                    # The "support" node is on the same column as the parent clade,
                    # and on the same row as the child clade. It is used to create the
                    # Edge config: parent -> support -> child

                    support_id = clade_id + 's' + str(n)
                    child_id = clade_id + 'c' + str(n)
                    pos_y_child = row_positions[child] * ylen

                    cy_support_node = {
                        'data': {'id': support_id},
                        'position': {'x': pos_x, 'y': pos_y_child},
                        'grabbable': grabbable,
                        'classes': 'support'
                    }

                    cy_support_edge = {
                        'data': {
                            'source': clade_id,
                            'target': support_id,
                            'sourceCladeId': clade_id
                        },
                    }

                    cy_edge = {
                        'data': {
                            'source': support_id,
                            'target': child_id,
                            'length': clade.branch_length,
                            'sourceCladeId': clade_id
                        },
                    }

                    if clade.confidence and clade.confidence.value:
                        cy_source['data']['confidence'] = clade.confidence.value

                    nodes.append(cy_support_node)
                    edges.extend([cy_support_edge, cy_edge])

                    add_to_elements(child, child_id)

            col_positions = get_col_positions(tree)
            row_positions = get_row_positions(tree)

            nodes = []
            edges = []

            add_to_elements(tree.clade, 'r')

            return nodes, edges

        tree = Phylo.read(tree_file, 'phyloxml')
        nodes, edges = generate_elements(tree)
        elements = nodes + edges

        layout = {'name': 'preset'}

        stylesheet = [
            {
                'selector': '.nonterminal',
                'style': {
                    'label': 'data(confidence)',
                    'background-opacity': 0,
                    "text-halign": "left",
                    "text-valign": "top",
                }
            },
            {
                'selector': '.support',
                'style': {'background-opacity': 0}
            },
            {
                'selector': 'edge',
                'style': {
                    "source-endpoint": "inside-to-node",
                    "target-endpoint": "inside-to-node",
                }
            },
            {
                'selector': '.terminal',
                'style': {
                    'label': 'data(name)',
                    'width': 10,
                    'height': 10,
                    "text-valign": "center",
                    "text-halign": "right",
                    'background-color': '#222222'
                }
            }
        ]

        app.layout = html.Div([
            cyto.Cytoscape(
                id='cytoscape-usage-phylogeny',
                elements=elements,
                stylesheet=stylesheet,
                layout=layout,
                style={
                    'height': '95vh',
                    'width': '100%'
                }
            )
        ])

        @app.callback(Output('cytoscape-usage-phylogeny', 'stylesheet'),
                      Input('cytoscape-usage-phylogeny', 'mouseoverEdgeData'))
        def color_children(edgeData):
            if edgeData is None:
                return stylesheet

            if 's' in edgeData['source']:
                val = edgeData['source'].split('s')[0]
            else:
                val = edgeData['source']

            children_style = [{
                'selector': 'edge[source *= "{}"]'.format(val),
                'style': {
                    'line-color': 'blue'
                }
            }]

            return stylesheet + children_style

        return super(PhyloTreeDetailView, self).get(request, *args, **kwargs)

    def get_context_data(self, **kwargs):
        context = super(PhyloTreeDetailView, self).get_context_data(**kwargs)
        context['bre'] = "Filogenetik Ağaç"
        return context


class AlignmentChartView(generic.DetailView):
    template_name = "bioinformatic/alignments/alignment_chart.html"
    model = MultipleSequenceAlignment

    def get(self, request, *args, **kwargs):
        app = DjangoDash('AlignmentChart')

        handle = os.path.join(BASE_DIR, "media", "MultipleSequenceAlignment", f"{self.request.user}",
                              f"{self.request.user}_out_alignment.fasta")

        data = open(handle).read()

        app.layout = html.Div([
            dashbio.AlignmentChart(
                id='alignment-viewer-eventDatum-usage',
                data=data,
            ),
            html.P('Yaptığınız İşlemler'),
            html.Div(id='alignment-viewer-eventDatum-usage-output')
        ])

        @app.callback(
            Output('alignment-viewer-eventDatum-usage-output', 'children'),
            Input('alignment-viewer-eventDatum-usage', 'eventDatum')
        )
        def update_output(value):
            if value is None:
                return 'Bir Sorun Oluştu'

            value = json.loads(value)

            if len(value.keys()) == 0:
                return 'Bir Sorun Oluştu'

            return [
                html.Div('- {}: {}'.format(key, value[key]))
                for key in value.keys()
            ]

        return super(AlignmentChartView, self).get(request, *args, **kwargs)

    def get_context_data(self, **kwargs):
        context = super(AlignmentChartView, self).get_context_data(**kwargs)
        context['bre'] = "Multiple Sekans Alignment Haritası"
        return context


class AlignmentClusterGramView(generic.DetailView):
    template_name = "bioinformatic/alignments/clustergram.html"
    model = MultipleSequenceAlignment

    def get(self, request, *args, **kwargs):
        if request.user.is_anonymous:
            from django.conf import settings
            messages.error(request, "Lütfen Giriş Yapınız")
            return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))

        return super(AlignmentClusterGramView, self).get(request, *args, **kwargs)

    def get_context_data(self, **kwargs):
        context = super(AlignmentClusterGramView, self).get_context_data(**kwargs)
        context['bre'] = "Multiple Sekans Alignment ClusterGram"
        handle = os.path.join(BASE_DIR, "media", "MultipleSequenceAlignment", f"{self.request.user}",
                              f"{self.request.user}_cluster.csv")

        df = pd.read_csv(handle).set_index('organizma')

        columns = list(df.columns.values)
        rows = list(df.index)

        clustergram = dashbio.Clustergram(
            data=df.loc[rows].values,
            row_labels=rows,
            column_labels=columns,

            color_threshold={
                'row': 250,
                'col': 700
            },
            height=800,
            width=800,
            color_list={
                'row': ['#636EFA', '#00CC96', '#19D3F3'],
                'col': ['#AB63FA', '#EF553B'],
                'bg': '#506784'
            },
            line_width=3
        )
        context['fig'] = clustergram.to_html()
        return context
