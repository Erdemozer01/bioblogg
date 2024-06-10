from django import forms


class ArticleForm(forms.Form):
    email = forms.EmailField()
    term = forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control'}))
