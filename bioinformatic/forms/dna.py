from django import forms


class DNASekansForm(forms.Form):
    dna = forms.CharField(label="DNA SEKANSI", widget=forms.Textarea())



