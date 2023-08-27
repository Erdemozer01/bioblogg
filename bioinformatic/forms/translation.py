from django import forms

TABLE = (
    ("", "------------"),
    ("1", "The Standard Code"),
    ("2", "The Vertebrate Mitochondrial Code"),
    ("3", "The Yeast Mitochondrial Code"),
    ("4", "The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code"),
    ("5", "The Invertebrate Mitochondrial Code"),
    ("6", "The Ciliate, Dasycladacean and Hexamita Nuclear Code"),
    ("9", "The Echinoderm and Flatworm Mitochondrial Code"),
    ("10", "The Euplotid Nuclear Code"),
    ("11", "The Bacterial, Archaeal and Plant Plastid Code, prokaryotic viruses"),
    ("12", "The Alternative Yeast Nuclear Code"),
    ("13", "The Ascidian Mitochondrial Code"),
    ("14", "The Alternative Flatworm Mitochondrial Code"),
    ("16", "Chlorophycean Mitochondrial Code"),
    ("21", "Trematode Mitochondrial Code"),
    ("22", "Scenedesmus obliquus Mitochondrial Code"),
    ("23", "Thraustochytrium Mitochondrial Code"),
    ("24", "Rhabdopleuridae Mitochondrial Code"),
    ("25", "Candidate Division SR1 and Gracilibacteria Code"),
    ("26", "Pachysolen tannophilus Nuclear Code"),
    ("27", "Karyorelict Nuclear Code"),
    ("28", "TCondylostoma Nuclear Code"),
    ("29", "Mesodinium Nuclear Code"),
    ("30", "Peritrich Nuclear Code"),
    ("31", "Blastocrithidia Nuclear Code"),
    ("33", "Cephalodiscidae Mitochondrial UAA-Tyr Code"),
)


class DNAFastaFileTranslateForm(forms.Form):
    file = forms.FileField(label="DNA Fasta Dosyası Giriniz")
    translate_table = forms.ChoiceField(choices=TABLE, label="Dönüşüm Tablosu Seçiniz")


class GenbankTranslationForm(forms.Form):
    table = forms.ChoiceField(
        label='Dönüşüm Tablosu',
        choices=TABLE,
        help_text='<p style="font-size: x-small"><b>Not: </b> Dönüşüm Tablosu verileri<a target="_blank" href="https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi"> NCBI </a> sitesinden alınmıştır. </p> ',
        widget=forms.Select(
            attrs={
                'class': 'form-control',
            }
        )
    )


class TranslationForm(forms.Form):
    table = forms.ChoiceField(
        label='Dönüşüm Tablosu',
        choices=TABLE,
        help_text=' <p style="font-size: x-small"><b>Not: </b> Dönüşüm Tablosu verileri \
        <a target="_blank" href="https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi"> \
        NCBI </a> sitesinden alınmıştır. </p> ',
        widget=forms.Select(
            attrs={
                'class': 'form-control',

            }
        )
    )
    sequence = forms.CharField(
        label='Sekans',
        widget=forms.Textarea(
            attrs={
                'class': 'form-control',
                'placeholder': 'Dna Sekansı Giriniz'

            }
        )
    )
