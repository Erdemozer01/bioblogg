from django import forms
from bioinformatic.models.bioinformatic import BioinformaticAnalizModel
from bioinformatic.choices import WRITE_FILE_FORMAT

class FastaWritingForm(forms.ModelForm):
    file_format = forms.ChoiceField(choices=WRITE_FILE_FORMAT, label="Dosya Formatı")

    class Meta:
        model = BioinformaticAnalizModel
        fields = ['molecule_id', 'description', 'seq', 'file_format']
