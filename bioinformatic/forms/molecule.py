from django import forms
from bioinformatic.models import MolecularModel


class MoleculeForm(forms.Form):
    file = forms.FileField(
        label="Dosya Seçiniz",
        help_text=".pdb veya .cif uzanlı olmalıdır. max. 2mb",
        required=False,
    )

    pdb_id = forms.CharField(max_length=10, min_length=4, widget=forms.TextInput(attrs={'placeholder': '1FAT'}),
                             label="PDB İD", help_text="Protein veritabanında id numarasına göre molekül görüntüleme yapabilirsiniz", required=False)


class MultipleMoleculeForm(forms.ModelForm):
    class Meta:
        model = MolecularModel
        fields = ['in_file']
        widgets = {
            'in_file': forms.MultipleHiddenInput(attrs={'multiple': True}),
        }

        help_texts = {
            'in_file': ".pdb veya .cif uzanlı olmalıdır. max. 2mb"
        }

        labels = {
            'in_file': "Dosya Seçiniz"
        }
