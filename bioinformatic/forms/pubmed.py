from django import forms
from bioinformatic.models import PubMedArticle


class PubMedForm(forms.ModelForm):
    class Meta:
        model = PubMedArticle
        fields = ['email', 'article_id']

        widgets = {
            'article_id': forms.TextInput(attrs={'placeholder': "PMID"})
        }

        labels = {
            'article_id': "Makale ID"
        }
