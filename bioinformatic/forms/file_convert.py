from django import forms
from bioinformatic.choices import *


class FileConvertForm(forms.Form):
    molecule = forms.ChoiceField(
        choices=MOLECULE_TYPE, label="Molekül tipi",
        widget=forms.Select(attrs={'class': 'custom-select'})
    )
    reading_file_format = forms.ChoiceField(
        choices=READ_FILE_FORMAT, label="Dönüştürülecek dosya formatı",
        widget=forms.Select(attrs={'class': 'custom-select'})
    )
    writing_file_format = forms.ChoiceField(
        choices=WRITE_FILE_FORMAT, label="Dönüştürmek istediğiniz dosya formatı",
        widget=forms.Select(attrs={'class': 'custom-select'})
    )

    file = forms.FileField(label="Dönüştürmek istediğiniz dosyayı seçiniz", widget=forms.FileInput())


