from django import forms
from bioinformatic.models import BioinformaticAnalizModel
from bioinformatic.choices import ALIGNMENT_MODE, MATRIS

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


class FileResulSelect(forms.ModelForm):
    class Meta:
        model = BioinformaticAnalizModel
        fields = ["select"]
        labels = {
            'select': ""
        }


class AlignResultForm(forms.Form):
    align = forms.Textarea()


class AlignmentForm(forms.Form):
    mode = forms.Select(choices=ALIGNMENT_MODE)
    matrix = forms.Select(choices=MATRIS)
    seq1 = forms.Textarea()
    seq2 = forms.Textarea()

