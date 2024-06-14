from django import forms
from bioinformatic.choices import ENTREZ_SELECT


class EntrezSelectForm(forms.Form):
    email = forms.EmailField()
    select = forms.ChoiceField(choices=ENTREZ_SELECT, widget=forms.Select(attrs={'class': 'form-control'}),
                               label="Yapmak istediğiniz işlem")


class ArticleForm(forms.Form):
    email = forms.EmailField()
    term = forms.CharField(
        widget=forms.TextInput(attrs={'class': 'form-control', 'placeholder': 'Aramak istediğiniz başlık yada terim'}),
        label="Terim")
