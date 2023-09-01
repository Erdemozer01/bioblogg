from django import forms

TABLE = (
    ("", "------------"),
    ("1", "Standart Kod"),
    ("2", "Omurgalı Mitokondri Kodu"),
    ("3", "Maya Mitokondri Kodu"),
    ("4", "Küf, Protozoon ve Kölenterat Mitokondri Kodu ve Mikoplazma/Spiroplasma Kodu"),
    ("5", "Omurgasız Mitokondri Kodu"),
    ("6", "Siliat, Dasikladas ve Heksamita Nükleer Kodu"),
    ("9", "Ekinoderm ve Yassı Solucan Mitokondri Kodu"),
    ("10", "Euplotid Nükleer Kodu"),
    ("11", "Bakteriyel, Arkeal ve Bitki Plastid Kodu, prokaryotik virüsler"),
    ("12", "Alternatif Maya Nükleer Kodu"),
    ("13", "Ascidian Mitokondri Kodu"),
    ("14", "Alternatif Yassı Solucan Mitokondri Kodu"),
    ("16", "Klorofis Mitokondri Kodu"),
    ("21", "Trematod Mitokondriyal Kodu"),
    ("22", "Scenedesmus obliquus Mitokondri Kodu"),
    ("23", "Thraustochytrium Mitokondri Kodu"),
    ("24", "Rhabdopleuridae Mitokondri Kodu"),
    ("25", "Aday Bölüm SR1 ve Gracilibacteria Kodu"),
    ("26", "Pachysolen tannophilus Nükleer Kodu"),
    ("27", "Karyorelict Nükleer Kodu"),
    ("28", "Kondilostoma Nükleer Kodu"),
    ("29", "Mezodinyum Nükleer Kodu"),
    ("30", "Peritrich Nükleer Kodu"),
    ("31", "Blastocrithidia Nükleer Kodu"),
    ("33", "Cephalodiscidae Mitokondriyal UAA-Tyr Kodu"),
)


class DNAFastaFileTranslateForm(forms.Form):

    translate_table = forms.ChoiceField(
        choices=TABLE,
        label="Dönüşüm Tablosu Seçiniz",
        widget=forms.Select(
            attrs={
                'class': 'custom-select',
            },
        ),

        required=True
    )

    file = forms.FileField(label="DNA Fasta Dosyası Giriniz", required=True)

    to_stop = forms.BooleanField(

        label="Stop kodonları dahil edilsin",

        required=False,

        widget=forms.CheckboxInput(attrs={
            'class': 'custom-control-input'
        })
    )


class GenbankTranslationForm(forms.Form):
    table = forms.ChoiceField(
        label='Dönüşüm Tablosu',
        choices=TABLE,
        help_text='<p style="font-size: x-small"><b>Not: </b> Dönüşüm Tablosu verileri<a target="_blank" href="https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi"> NCBI </a> sitesinden alınmıştır. </p> ',
        widget=forms.Select(
            attrs={
                'class': 'custom-select',
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
                'class': 'custom-select',

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

    to_stop = forms.BooleanField(

        label="Stop kodonları dahil edilsin",

        required=False,

        widget=forms.CheckboxInput(attrs={
            'class': 'custom-control-input'
        })
    )
