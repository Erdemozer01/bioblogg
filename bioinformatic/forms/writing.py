from django import forms


class FastaWritingForm(forms.Form):
    id = forms.CharField(
        label="id",
        widget=forms.TextInput(
            attrs={
                'class': 'form-control',
                'placeholder': 'id numarası'
            }
        )
    )

    description = forms.CharField(
        label="Description",
        widget=forms.TextInput(
            attrs={
                'class': 'form-control',
                'placeholder': 'Tanım'
            }
        )
    )

    sequence = forms.CharField(
        label="Sequence",
        widget=forms.Textarea(
            attrs={
                'class': 'form-control',
                'placeholder': 'Sekans'
            }
        )
    )

class GenbankWritingForm(forms.Form):
    lokus = forms.CharField(
        required=True,
        widget=forms.TextInput(
            attrs={
                'class': 'form-control',
                'placeholder': 'Sequence name, e.g. gene name, LOKUS'
            }
        )
    )

    molecule_type = forms.CharField(
        required=True,
        widget=forms.TextInput(
            attrs={
                'class': 'form-control',
                'placeholder': 'Molekül tipi. Ex. DNA'
            }
        )
    )

    keywords = forms.CharField(
        required=False,
        widget=forms.TextInput(

            attrs={
                'class': 'form-control',
                'placeholder': 'Anahtar Kelimeler'
            }
        )
    )

    taxonomy = forms.CharField(
        widget=forms.TextInput(
            attrs={
                'class': 'form-control',
                'placeholder': 'Taksonomi'
            }
        )
    )

    references = forms.CharField(
        widget=forms.TextInput(
            attrs={
                'class': 'form-control',
                'placeholder': 'Referans'
            }
        )
    )

    sequence = forms.CharField(
        widget=forms.Textarea(
            attrs={
                'class': 'form-control',
                'placeholder': 'Sekans'
            }
        )
    )

    id = forms.CharField(
        widget=forms.TextInput(
            attrs={
                'class': 'form-control',
                'placeholder': 'id numarası'
            }
        )
    )

    description = forms.CharField(
        widget=forms.TextInput(
            attrs={
                'class': 'form-control',
                'placeholder': 'Tanımlama'
            }
        )
    )

    dbxref = forms.CharField(
        required=False,
        help_text='<small style="font-size: small; font-style: italic; '
                  'margin-left: 10px">Veri tabanı referansı</small>',
        widget=forms.TextInput(
            attrs={
                'class': 'form-control',
                'placeholder': 'dbxref'
            }
        )
    )

    source = forms.CharField(
        required=True,
        widget=forms.TextInput(
            attrs={
                'class': 'form-control',
                'placeholder': 'source'
            }
        )
    )

    organism = forms.CharField(
        required=True,
        widget=forms.TextInput(
            attrs={
                'class': 'form-control',
                'placeholder': 'Organizma'
            }
        )
    )

    features = forms.CharField(
        required=True,
        widget=forms.TextInput(
            attrs={
                'class': 'form-control',
                'placeholder': 'Özellik'
            }
        )
    )

