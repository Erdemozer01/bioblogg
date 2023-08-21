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


class ContactProfileForm(forms.ModelForm):
    class Meta:
        model = ContactModel
        fields = ['title', 'contact_email', 'content']

        widgets = {
            'content': forms.Textarea(
                attrs={'class': "form-control border border-light", 'placeholder': 'Mesajınız'},
            ),

            'contact_email': forms.EmailInput(attrs={'class': 'form-control', 'placeholder': 'Email adresiniz'}),
            'title': forms.TextInput(attrs={'class': 'form-control', 'placeholder': 'Konu'}),
        }

        labels = {
            'contact_email': "Email",
            'title': "Konu",
        }
