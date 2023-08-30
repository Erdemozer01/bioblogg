from django import forms
from bioinformatic.models import FastaRead, GenbankRead, FileUploadModel

FILE_TYPE = (
    ("fasta", "fasta"),
    ("genbank", "genbank"),
)

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


# creating a form

class FileTypeSelect(forms.Form):
    geeks_field = forms.ChoiceField(choices=FILE_TYPE)


class FileReadForm(forms.Form):
    org_type = forms.ChoiceField(choices=TABLE, label="Dönüşüm Tablosu",
                                 widget=forms.Select(attrs={'class': 'custom-select'}))
    file = forms.FileField(
        label='Dosya Seçiniz',
        help_text="Not : Max. dosya boyutu 5mb olmalıdır.",
        widget=forms.FileInput()
    )


class FastaIdForm(forms.Form):
    gene = forms.ModelChoiceField(queryset=FastaRead.objects.all(), label='Gen Bölgesi Seçiniz')


class GenbankIdForm(forms.Form):
    gene = forms.ModelChoiceField(queryset=GenbankRead.objects.all(), label='Gen Bölgesi Seçiniz')


class XmlIdForm(forms.Form):
    hit_id = forms.Textarea()


class FileUploadModelForm(forms.ModelForm):
    class Meta:
        model = FileUploadModel
        fields = ['file']


class MultipleUploadFileForm(forms.Form):
    file_field = forms.FileField(widget=forms.MultipleHiddenInput(
        attrs={'multiple': True}),
        label="Fasta Dosyaları")
