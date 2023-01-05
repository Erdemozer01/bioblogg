from django import forms
from accounts.models import ContactModel


class ContactForm(forms.ModelForm):
    accept = forms.BooleanField()

    class Meta:
        model = ContactModel
        fields = ['content', 'contact_email']

        widgets = {
            'content': forms.Textarea(
                attrs={'class': "form-control border border-light", 'placeholder': 'Mesajınız'},
            ),

            'contact_email': forms.EmailInput(attrs={'class': 'form-control', 'placeholder': 'Email adresiniz'})
        }
