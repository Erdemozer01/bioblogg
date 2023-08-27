from django import forms
from bioinformatic.models import FastaRead, GenbankRead, FileUploadModel

FILE_TYPE = (
    ("fasta", "fasta"),
    ("genbank", "genbank"),
)


# creating a form

class FileTypeSelect(forms.Form):
    geeks_field = forms.ChoiceField(choices=FILE_TYPE)


class FileReadForm(forms.Form):
    file = forms.FileField(label='Dosya Seçiniz')


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
