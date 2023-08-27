from django.shortcuts import render
from bioinformatic.forms.translation import TranslationForm
import Bio.Data.CodonTable
from Bio.Seq import Seq
from Bio.SeqUtils import GC


def translation(request):
    form = TranslationForm(request.POST or None)
    if request.method == "POST":
        if form.is_valid():
            table = form.cleaned_data['table']
            sequence = form.cleaned_data['sequence'].upper()
            sequence = sequence.replace(" ", "")
            sequence = sequence.replace("\n", "")
            sequence = sequence.replace("\t", "")
            sequence = sequence.replace("\r", "")
            sequence = sequence.replace("0", "")
            sequence = sequence.replace("1", "")
            sequence = sequence.replace("2", "")
            sequence = sequence.replace("3", "")
            sequence = sequence.replace("4", "")
            sequence = sequence.replace("5", "")
            sequence = sequence.replace("6", "")
            sequence = sequence.replace("7", "")
            sequence = sequence.replace("8", "")
            sequence = sequence.replace("9", "")

            complement = Seq("{}".format(sequence)).complement()
            revese = Seq("{}".format(sequence)).reverse_complement()
            transcribe = Seq("{}".format(sequence)).transcribe()

            try:

                translate = Seq("{}".format(sequence)).translate(table='{}'.format(table))

            except Bio.Data.CodonTable.TranslationError:

                msg = "Codon Tablo Hatası"

                return render(request, 'accounts/404.html', {'msg': msg})

            return render(request, 'bioinformatic/translation/result.html',
                          context={
                              'sequence': sequence,
                              'len': len(sequence),
                              'GC': GC(sequence),
                              'table': table,
                              'complement': complement,
                              'revese': revese,
                              'transcribe': transcribe,
                              'translate': translate,
                              'bre': 'Protein Sentezi'
                          })
        else:
            error = 'hata'
            return render(request, 'accounts/404.html', context={
                'msg': "HATA"
            })

    return render(request, 'bioinformatic/translation/analize.html',
                  context={
                      'form': form,
                      "bre": "Protein Sentezi"
                  })
