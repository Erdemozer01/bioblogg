from django import forms


class DNASekansForm(forms.Form):
    dna = forms.CharField(label="", widget=forms.Textarea(attrs={'placeholder': 'DNA SEKANSI GİRİNİZ'}))


class SequenceSlicingForm(forms.Form):
    seq = forms.CharField(label="SEKANS", widget=forms.Textarea(attrs={'placeholder': 'SEKANS GİRİNİZ'}))
    start = forms.IntegerField(label="Başlangıç", min_value=1)
    finish = forms.IntegerField(label="Bitiş", min_value=1)
