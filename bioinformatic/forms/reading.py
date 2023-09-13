from django import forms
from bioinformatic.models import BioinformaticAnalizModel
from bioinformatic.choices import ALIGNMENT_MODE, MATRIS
from ckeditor.widgets import CKEditorWidget

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



class AlignmentForm(forms.Form):
    mode = forms.ChoiceField(choices=ALIGNMENT_MODE, widget=forms.Select(attrs={
        'class': 'custom-select'
    }))
    matrix = forms.ChoiceField(choices=MATRIS, widget=forms.Select(attrs={
        'class': 'custom-select'
    }))
    seq1 = forms.CharField(widget=forms.Textarea(attrs={'placeholder': 'Sekans 1'}), label="Sekans 1")
    seq2 = forms.CharField(widget=forms.Textarea(attrs={'placeholder': 'Sekans 2'}), label="Sekans 2")
