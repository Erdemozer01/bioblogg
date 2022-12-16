from django import forms
from blog.models import Posts, Category


class CreatePostModelForm(forms.ModelForm):
    class Meta:
        model = Posts
        exclude = ('author',)

        widgets = {
            'cover': forms.ClearableFileInput(attrs={'class': 'form-control', 'type': 'file'}),
            'extra_image1': forms.ClearableFileInput(attrs={'class': 'form-control', 'type': 'file'}),
            'extra_image2': forms.ClearableFileInput(attrs={'class': 'form-control', 'type': 'file'}),
            'extra_image3': forms.ClearableFileInput(attrs={'class': 'form-control', 'type': 'file'}),
            'title': forms.TextInput(attrs={'class': 'form-control', 'placeholder': 'Gönderi Başlığı'}),
            'category': forms.Select(attrs={'class': 'form-control'}),
            'status': forms.Select(attrs={'class': 'form-control'}),
        }

        labels = {
            'cover': "Kapak Fotosu",
            'category': "Kategori Seçiniz",
            'title': "Gönderi Başlığı",
            'text': ""
        }


class CreateCategoryModelForm(forms.ModelForm):
    class Meta:
        model = Category
        fields = ['image', 'title', 'explain']

        widgets = {
            'image': forms.ClearableFileInput(attrs={'class': 'custom-file-input', 'type': 'file'}),
            'title': forms.TextInput(attrs={'class': 'form-control', 'placeholder': 'Kategori Adı'}),
            'explain': forms.Textarea(attrs={'class': 'form-control', 'placeholder': 'Kategori Açıklaması'}),
        }
