from django import forms
from bioinformatic.models import BioinformaticModel, FileModel
from bioinformatic.choices import *


class FileReadingForm(forms.ModelForm):
    file = forms.FileField(
        label='Dosya',
        help_text="Not : Max. dosya boyutu 30 mb olmalıdır.",
    )


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
                    'class': 'form-control'
                }
            ),
            'file': forms.FileInput(attrs={'class': 'custom-file-input'})
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


class MultipleSeqAlignmentFileForm(forms.Form):
    alignment_tools = forms.ChoiceField(label="Alingmen Araçları", choices=METHOD, widget=forms.Select(attrs={
        'class': 'custom-select'
    }))
    tree_alg = forms.ChoiceField(label="Algoritma", choices=TREE_ALGORITMA, widget=forms.Select(attrs={
        'class': 'custom-select'
    }))
    file = forms.FileField(label='Fasta Dosyası', help_text="Not : Max. dosya boyutu 25mb olmalıdır.")


class BlastForm(forms.Form):
    type = forms.ChoiceField(choices=RECORD_FORMAT, label="Blast Türü Seçiniz", widget=forms.Select(attrs={
        'class': 'custom-select'
    }))
    file = forms.FileField(label="Fasta Dosyası", required=False)
    gi = forms.CharField(label="Gİ (GENİNFO) Numarası", required=False)
    program = forms.ChoiceField(choices=BLAST_PROGRAM, label="Blast Programı Seçiniz", widget=forms.Select(attrs={
        'class': 'custom-select'
    }))
    database = forms.ChoiceField(choices=BLAST_DATABASE, label="Blast Veritabanı Seçiniz", widget=forms.Select(attrs={
        'class': 'custom-select'
    }))


class MultipleFileInput(forms.ClearableFileInput):
    allow_multiple_selected = True


class MultipleFileField(forms.FileField):
    def __init__(self, *args, **kwargs):
        kwargs.setdefault("widget", MultipleFileInput())
        super().__init__(*args, **kwargs)

    def clean(self, data, initial=None):
        single_file_clean = super().clean
        if isinstance(data, (list, tuple)):
            result = [single_file_clean(d, initial) for d in data]
        else:
            result = [single_file_clean(data, initial)]
        return result


class MultiMoleculeViewForm(forms.Form):
    files = MultipleFileField(
        label='Dosya yada dosyaları seçiniz',
        help_text="Not : Dosyanız pdb uzanlı olmalıdır. Dosya boyutu en fazla 15mb olmalıdır.",
    )


class SingleMoleculeViewForm(forms.Form):
    file = forms.FileField(
        label='DOSYA SEÇİNİZ',
        help_text="Not : Dosyanız pdb yada cif uzanlı olmalıdır. Dosya boyutu en fazla 9mb olmalıdır.",
    )
