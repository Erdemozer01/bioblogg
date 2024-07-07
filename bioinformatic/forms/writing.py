from datetime import datetime

from django import forms
from bioinformatic.models.bioinformatic import BioinformaticModel
from bioinformatic.choices import MOLECULE_TYPE, TOPOLOGY


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


class FastaWritingForm(forms.Form):
    molecule = forms.ChoiceField(choices=MOLECULE_TYPE, label="Molekül*", required=True,
                                 widget=forms.Select(attrs={'class': 'custom-select', 'placeholder': 'Molekül'}))
    molecule_id = forms.CharField(required=True, label="Molekül İd*",
                                  widget=forms.TextInput(attrs={'placeholder': 'Molekül İD'}))
    name = forms.CharField(required=True, label="Ad*",
                           widget=forms.TextInput(attrs={'placeholder': 'Molekül Adı'}))
    description = forms.CharField(label="Tanım*", required=True,
                                  widget=forms.TextInput(attrs={'placeholder': 'Tanımlama'}))
    sequence = forms.CharField(label="Sekans*", required=True,
                               widget=forms.Textarea(attrs={'placeholder': 'Sekans giriniz'}))


class GenbankWritingForm(forms.Form):
    molecule = forms.ChoiceField(choices=MOLECULE_TYPE, label="Molekül*", required=True,
                                 widget=forms.Select(attrs={'class': 'custom-select', 'placeholder': 'Molekül'}))
    molecule_id = forms.CharField(required=True, label="Molekül İd*",
                                  widget=forms.TextInput(attrs={'placeholder': 'Molekül İd'}))
    name = forms.CharField(required=True, label="Ad*".title(),
                           widget=forms.TextInput(attrs={'placeholder': 'Ad'}))
    description = forms.CharField(label="Tanım*", required=True,
                                  widget=forms.TextInput(attrs={'placeholder': 'Tanım'}))
    db_xrefs = forms.CharField(label="Veri tabanı referans*", required=True,
                               widget=forms.TextInput(attrs={'placeholder': 'Veri tabanı referans'}))

    # Anatasyonlar
    taxonomy = forms.CharField(label="Taksonomi*", widget=forms.TextInput(attrs={'placeholder': 'Taksonomi'}))
    source = forms.CharField(label="Kaynak*", widget=forms.TextInput(attrs={'placeholder': 'Kaynak'}))
    organism = forms.CharField(label="Organizma*", widget=forms.TextInput(attrs={'placeholder': 'Organizma'}))
    accessions = forms.CharField(label="Erişim Numarası*",
                                 widget=forms.TextInput(attrs={'placeholder': 'Erişim Numarası'}))
    topology = forms.ChoiceField(label="Topoloji*", choices=TOPOLOGY,
                                 widget=forms.Select(attrs={'class': 'custom-select'}))
    # Sekans
    sequence = forms.CharField(label="Sekans*", required=True,
                               widget=forms.Textarea(attrs={'placeholder': 'Sekans giriniz'.title()}))
