from django import forms
from django.contrib.auth.models import User
from django.contrib.auth.forms import UserCreationForm
from .models import Profile


class UserRegistrationForm(UserCreationForm):
    terms = forms.BooleanField(required=True, label="Kullanıcı Şartlarını kabul ediyorum.")

    class Meta:
        model = User
        fields = ['username', 'first_name', 'last_name', 'email', 'password1', 'password2', 'terms']

    def clean_email(self):
        data = self.cleaned_data['email']
        if User.objects.filter(email=data).exists():
            raise forms.ValidationError('Bu Email zaten kullanılmakta')
        return data


class UserEditForm(forms.ModelForm):
    class Meta:
        model = User
        fields = ['username', 'first_name', 'last_name', 'email', 'is_superuser', 'is_staff', 'is_active']
        widgets = {
            'username': forms.TextInput(attrs={'class': 'form-control', 'placeholder': 'Kullanıcı Adı'}),
            'first_name': forms.TextInput(attrs={'class': 'form-control', 'placeholder': 'Adı'}),
            'last_name': forms.TextInput(attrs={'class': 'form-control', 'placeholder': 'Soyadı:'}),
            'email': forms.EmailInput(attrs={'class': 'form-control', 'placeholder': 'Email:'}),
        }


class UserProfileEditForm(forms.ModelForm):
    class Meta:
        model = User
        fields = ['username', 'first_name', 'last_name', 'email']
        widgets = {
            'username': forms.TextInput(attrs={'class': 'form-control', 'placeholder': 'Kullanıcı Adı'}),
            'first_name': forms.TextInput(attrs={'class': 'form-control', 'placeholder': 'Adı'}),
            'last_name': forms.TextInput(attrs={'class': 'form-control', 'placeholder': 'Soyadı'}),
            'email': forms.EmailInput(attrs={'class': 'form-control', 'placeholder': 'Email'}),
        }


class ProfileEditForm(forms.ModelForm):
    class Meta:
        model = Profile
        fields = ['cover', 'avatar', 'phone', 'location', 'job', 'about', 'facebook', 'twitter', 'instagram']

