from django import forms
from bioinformatic.models import FileFormat



molecule_type = (
    ('', '-' * 30),
    ('DNA', 'DNA'),
    ('RNA', 'RNA'),
    ('protein', 'PROTEİN'),
)


class FileConvertForm(forms.Form):
    in_file = forms.FileField(label="Dosya")
    molecule_type = forms.ChoiceField(label="Molekül Tipi", choices=molecule_type)



class FileConvertModelForm(forms.ModelForm):
    class Meta:
        model = FileFormat
        fields = ['name']
