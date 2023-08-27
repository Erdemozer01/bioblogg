from django import forms
from django.conf import settings

METHOD = (
    ('', '------------'),
    ('nj', 'Neighbor Joining'),
    ('upgma', 'UPGMA'),
)


class PhyloGeneticTreeForm(forms.Form):
    method = forms.ChoiceField(choices=METHOD, label="Algoritme Seçiniz")
    files = forms.FileField(label="Fasta Dosyası Seçiniz", widget=forms.ClearableFileInput())

