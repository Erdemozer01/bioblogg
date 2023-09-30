from django import forms
from bioinformatic.models import BioinformaticModel
from bioinformatic.choices import ALIGNMENT_MODE, MATRIS


class FileReadingForm(forms.ModelForm):
    file = forms.FileField(label='Dosya', help_text="Not : Max. dosya boyutu 25 mb olmalıdır.")

    class Meta:
        model = BioinformaticModel
        fields = ["reading_file_format", 'molecule', "file"]
        widgets = {
            'reading_file_format': forms.Select(
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
        model = BioinformaticModel
        fields = ["trans_table", 'to_stop']
        widgets = {
            'trans_table': forms.Select(
                attrs={
                    'class': 'custom-select'
                })
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
