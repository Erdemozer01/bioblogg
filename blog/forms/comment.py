from django import forms
from blog.models import Comments
from ckeditor.widgets import CKEditorWidget


class AddCommentForm(forms.ModelForm):
    class Meta:
        model = Comments
        fields = ['comment']
        widgets = {
            'comment': CKEditorWidget(
                attrs={'class': "form-control border border-light",
                       "placeholder": "Yorumunuzu Yazın",
                       "data-bind-characters-target": "#charactersRemaining",
                       "label": "Yorum Yap",
                       "id": "CommentFormControlTextarea"
                       },
                config_name="blog"
            )
        }
