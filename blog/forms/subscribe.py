from django import forms
from blog.models import Subscribe


class SubscribeModelForm(forms.ModelForm):
    class Meta:
        model = Subscribe
        fields = "__all__"
