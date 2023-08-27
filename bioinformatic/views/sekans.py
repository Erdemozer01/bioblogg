from django.shortcuts import render
from Bio.Seq import Seq
from bioinformatic.forms.dna import DNASekansForm
from Bio.SeqUtils import GC


def sekans(request):
    form = DNASekansForm(request.POST or None)
    if request.method == "POST":
        if form.is_valid():

            sekans = form.cleaned_data['dna'].upper()

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

            my_seq = Seq("{}".format(sekans))

            protein = ['M', 'R', 'N', '	D', 'E', 'Q', 'H']

            istenmeyen = ['\n', '\r', '\t', '1', '2', '3', '4', '5']

            for aa in protein:
                for dna in my_seq:
                    if aa in dna:
                        return render(request, 'accounts/404.html', {'msg': 'Sekensda protein tespit edilmiştir.'})

            return render(request, 'bioinformatic/sekans/result.html', {'len': len(my_seq), 'A': my_seq.count('A'), 'G': my_seq.count('G'), 'C': my_seq.count('C'),
                           'T': my_seq.count('T'), 'GC': GC(my_seq), 'sekans': sekans})

        else:
            form = DNASekansForm(request.POST)

    return render(request, 'bioinformatic/sekans/analize.html', {'form': form})
