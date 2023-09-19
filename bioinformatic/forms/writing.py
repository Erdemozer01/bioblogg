from django import forms
from bioinformatic.models.bioinformatic import BioinformaticAnalizModel


class SelectWritingForm(forms.ModelForm):
    class Meta:
        model = BioinformaticAnalizModel
        fields = ['writing_file_format']
        widgets = {
            'writing_file_format': forms.Select(
                attrs={
                    'class': 'custom-select'
                }),
        }


class FileWritingForm(forms.ModelForm):
    molecule_id = forms.CharField(required=True, label="İD")

    class Meta:
        model = BioinformaticAnalizModel
        fields = [
            'molecule',
            'name',
            'molecule_id',
            'description',
            'dbxrefs',
            'source',
            'keywords',
            'organism',
            'accession',
            'seq'
        ]
        widgets = {
            'molecule': forms.Select(
                attrs={
                    'class': 'custom-select'
                }
            )
        }


class ConvertFile(forms.ModelForm):
    class Meta:
        model = BioinformaticAnalizModel
        fields = []