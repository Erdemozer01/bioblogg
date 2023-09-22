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


class FileWritingForm(forms.Form):
    molecule = forms.ChoiceField(choices=MOLECULE_TYPE, label="Molekül*", required=True,
                                 widget=forms.Select(attrs={'class': 'custom-select', 'placeholder': 'Molekül'}))
    molecule_id = forms.CharField(required=True, label="Molekül İD*".title(),
                                  widget=forms.TextInput(attrs={'placeholder': 'Molekül İD'}))
    name = forms.CharField(required=True, label="name*".title(),
                           widget=forms.TextInput(attrs={'placeholder': 'Molekül Adı'}))
    description = forms.CharField(label="description*".title(), required=True,
                                  widget=forms.TextInput(attrs={'placeholder': 'Tanımlama'}))
    db_xrefs = forms.CharField(label="db_xrefs".title())
    annotations = forms.CharField(label="annotations".title(), help_text=" Genbank dosyası için gerekli.")
    source = forms.CharField(label="source".title(), help_text="Genbank dosyası için gerekli.")
    organism = forms.CharField(label="organism".title(), help_text="Genbank dosyası için gerekli.")
    keywords = forms.CharField(label="keywords".title(), help_text="Genbank dosyası için gerekli.")
    accessions = forms.CharField(label="accessions".title(), help_text="Genbank dosyası için gerekli.")
    sequence = forms.CharField(label="Sequence*".title(), required=True,
                               widget=forms.Textarea(attrs={'placeholder': 'Sekans giriniz'.title()}))


class ConvertFile(forms.ModelForm):
    class Meta:
        model = BioinformaticModel
        fields = []
