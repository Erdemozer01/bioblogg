from django import forms
from django.forms.models import inlineformset_factory, modelform_factory, modelformset_factory

DIAGRAM_SHAPE = (
    ('', '------------'),
    ('ARROW', 'OK'),
    ('BOX', 'KUTU'),
)

DIAGRAM_FORMAT = (
    ('', '------------'),
    ('linear', 'DÜZ'),
    ('circular', 'DAİRESEL'),
)

from bioinformatic.models import RestrictionModel, DiagramModel


class GenomeDiagramForm(forms.Form):
    file = forms.FileField(label="Genbank Dosyası", help_text="Çoklu kayıt girmeyiniz")
    diagram_shape = forms.ChoiceField(label="Diagram Şekli", choices=DIAGRAM_SHAPE)
    diagram_format = forms.ChoiceField(label="Diagram Formatı", choices=DIAGRAM_FORMAT)
    fragment = forms.IntegerField(label="Fragment Sayısı", initial=1, min_value=1)


class GenomeDiagramAddEnzymesForm(forms.Form):
    enzymes = forms.CharField(label="Enzim Adı")
    site = forms.CharField(label="Bölge")


class RestrictionModelForms(forms.ModelForm):
    class Meta:
        model = DiagramModel
        exclude = ('user', 'out_file')


RestrictionModelFormFactory = modelformset_factory(RestrictionModel, exclude=('user',), extra=0, can_delete=True)
RestrictionModelFormFactory2 = modelformset_factory(RestrictionModel, exclude=('user',), extra=0)

RestrictionModelFormset = inlineformset_factory(
    DiagramModel, RestrictionModel, exclude=('user',)
)
