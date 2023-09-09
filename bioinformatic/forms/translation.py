from django import forms
from bioinformatic.choices import TRANS_TABLE


class TranslationForm(forms.Form):
    table = forms.ChoiceField(
        label='Dönüşüm Tablosu',
        choices=TRANS_TABLE,
        help_text=' <p style="font-size: x-small"><b>Not: </b> Dönüşüm Tablosu verileri \
        <a target="_blank" href="https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi"> \
        NCBI </a> sitesinden alınmıştır. </p> ',
        widget=forms.Select(
            attrs={
                'class': 'custom-select',

            }
        )
    )
    sequence = forms.CharField(
        label='Sekans',
        widget=forms.Textarea(
            attrs={
                'class': 'form-control',
                'placeholder': 'Dna Sekansı Giriniz'

            }
        )
    )

    to_stop = forms.BooleanField(

        label="Stop kodonları dahil edilsin",

        required=False,

        widget=forms.CheckboxInput(attrs={
            'class': 'custom-control-input'
        })
    )
