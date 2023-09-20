from django import forms
from bioinformatic.models.bioinformatic import BioinformaticModel, RecordModel
from bioinformatic.choices import MOLECULE_TYPE


class SelectWritingForm(forms.ModelForm):
    class Meta:
        model = BioinformaticModel
        fields = ['writing_file_format']
        widgets = {
            'writing_file_format': forms.Select(
                attrs={
                    'class': 'custom-select'
                }),
        }


class FileWritingForm(forms.ModelForm):
    molecule = forms.ChoiceField(choices=MOLECULE_TYPE, label="Molekül",
                                 widget=forms.Select(attrs={'class': 'custom-select'}))
    molecule_id = forms.CharField(required=True, label="İD")
    name = forms.CharField(required=True, label="name")
    db_xrefs = forms.CharField(label="db_xrefs")
    annotations = forms.CharField(label="annotations")
    source = forms.CharField(label="source")
    keywords = forms.CharField(label="keywords")
    accessions = forms.CharField(label="accessions")
    sequence = forms.CharField(label="sequence",
                               widget=forms.Textarea(attrs={'placeholder': 'Sekans giriniz'}))


class ConvertFile(forms.ModelForm):
    class Meta:
        model = BioinformaticModel
        fields = []
