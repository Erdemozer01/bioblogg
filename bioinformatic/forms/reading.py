from django import forms
from bioinformatic.models import BioinformaticAnalizModel


class FileReadingForm(forms.ModelForm):
    class Meta:
        model = BioinformaticAnalizModel
        fields = ["file_format", 'molecule', "read_file"]
        widgets = {
            'file_format': forms.Select(
                attrs={
                    'class': 'custom-select'
                }),
            'molecule': forms.Select(
                attrs={
                    'class': 'custom-select'
                }
            )
        }


class TranslateForm(forms.ModelForm):
    class Meta:
        model = BioinformaticAnalizModel
        fields = ["trans_table", 'to_stop']
        widgets = {
            'trans_table': forms.Select(
                attrs={
                    'class': 'custom-select'
                })
        }
