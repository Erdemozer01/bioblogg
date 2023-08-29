from django import forms


class DNASekansForm(forms.Form):
    dna = forms.CharField(label="", widget=forms.Textarea(attrs={'placeholder': 'DNA SEKANSI GİRİNİZ'}))
