from django import forms


class AddFastaData(forms.Form):
    file = forms.FileField(label='Dosya Seçiniz')

    fasta_id = forms.CharField(
        label="İD",
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