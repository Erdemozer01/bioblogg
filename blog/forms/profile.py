from django import forms
from accounts.models import Profile,SocialMedia


class ProfileModelForm(forms.ModelForm):
    class Meta:
        model = Profile
        exclude = ('user', 'cover')
        widgets = {
            'birth_day': forms.SelectDateWidget(attrs={'class': 'form-control'}, years=range(1900, 2050)),
            'first_name': forms.TextInput(attrs={'class': 'form-control'}),
            'avatar': forms.ClearableFileInput(attrs={'class': 'form-control'}),
            'gender': forms.Select(attrs={'class': 'form-control'}),
        }


class SocialMediaModelForm(forms.ModelForm):
    class Meta:
        model = SocialMedia
        exclude = ('user',)
