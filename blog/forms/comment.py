from django import forms
from blog.models import Comments


class AddCommentForm(forms.ModelForm):
    class Meta:
        model = Comments
        fields = ['comment']
        widgets = {
            'comment': forms.Textarea(
                attrs={'class': "form-control border border-light",
                       "placeholder": "Yorumunuzu Yazın",
                       "data-bind-characters-target": "#charactersRemaining",
                       "label": "Yorum Yap",
                       "id": "CommentFormControlTextarea"
                       },
            )
        }
