from django import forms
from bioinformatic.models import MedlineArticle


class MedlineArticleForm(forms.ModelForm):
    class Meta:
        model = MedlineArticle
        fields = ['email', 'term']

        widgets = {
            'term': forms.TextInput(attrs={'placeholder': "Aranacak Kelime"})
        }

        labels = {
            'term': "Aranacak Kelime"
        }
