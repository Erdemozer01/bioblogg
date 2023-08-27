from django import forms
from bioinformatic.models import MultipleSequenceAlignment

MATRIS = (
    ('', '------------'),
    ('BENNER22', 'BENNER22'),
    ('BENNER6', 'BENNER6'),
    ('BENNER74', 'BENNER74'),
    ('BLOSUM45', 'BLOSUM45'),
    ('BLOSUM50', 'BLOSUM50'),
    ('BLOSUM62', 'BLOSUM62'),
    ('BLOSUM80', 'BLOSUM80'),
    ('BLOSUM90', 'BLOSUM90'),
    ('DAYHOFF', 'DAYHOFF'),
    ('FENG', 'FENG'),
    ('HOXD70', 'HOXD70'),
    ('JOHNSON', 'JOHNSON'),
    ('JONES', 'JONES'),
    ('LEVIN', 'LEVIN'),
    ('MCLACHLAN', 'MCLACHLAN'),
    ('MDM78', 'MDM78'),
    ('PAM250', 'PAM250'),
    ('PAM30', 'PAM30'),
    ('PAM70', 'PAM70'),
    ('RAO', 'RAO'),
    ('RISLER', 'RISLER'),
    ('SCHNEIDER', 'SCHNEIDER'),
    ('STR', 'STR'),
    ('TRANS', 'TRANS'),
)

METHOD = (
    ('', '------------'),
    ('muscle', 'MUSCLE'),
    ('clustalw2', 'CLUSTALW2'),
    ('omega', 'ClustalOmega'),
    ('paml', 'Maximum Likelihood (PAML)'),
)

ALGORITMA = (
    ('', '------------'),
    ('nj', 'Neighbor Joining'),
    ('upgma', 'UPGMA'),
)

MOLECULE_TYPE = (
    ('', '------------'),
    ('DNA', 'DNA'),
    ('RNA', 'RNA'),
    ('protein', 'PROTEİN'),
)

ALIGNMENT_FILE_TYPE = (
    ('', '------------'),
    ('fasta', 'FASTA'),
    ('clustal', 'CLUSTAL'),
    ('phylip', 'PHYLİB'),
    ('nexus', 'NEXUS'),
)

PALM_TOOLS = (
    ('', '------------'),
    ('baseml', 'BASEML'),
    ('basemlg', 'BASEMLG'),
    ('codeml', 'CODEML'),
)


class LocalForm(forms.Form):
    algo = forms.ChoiceField(
        label='Dizileme Matrisi Seçiniz',
        choices=MATRIS,
        widget=forms.Select(
            attrs={
                'class': 'form-control',

            }
        )
    )

    seq1 = forms.CharField(
        label="Sekans1 Giriniz",
        widget=forms.Textarea(
            attrs={
                'class': 'form-control',
                'placeholder': 'Sequence1'
            }
        )
    )

    seq2 = forms.CharField(
        label="Sekans2 Giriniz",
        widget=forms.Textarea(
            attrs={
                'class': 'form-control',
                'placeholder': 'Sequence2'

            }
        )
    )


class GlobalForm(forms.Form):
    seq1 = forms.CharField(
        label="Sekans1 Giriniz",
        widget=forms.Textarea(
            attrs={
                'class': 'form-control',
                'placeholder': 'Sequence1'
            }
        )
    )

    seq2 = forms.CharField(
        label="Sekans2 Giriniz",
        widget=forms.Textarea(
            attrs={
                'class': 'form-control',
                'placeholder': 'Sequence2'

            }
        )
    )


class MultipleSequenceAlignmentSelectForm(forms.Form):
    method = forms.ChoiceField(choices=METHOD, label="Multiple Sekans Alignment Aracı")


class MultipleSequenceAlignmentForm(forms.Form):
    file = forms.FileField(label="Fasta Dosyası Seçiniz")
    molecule_type = forms.ChoiceField(choices=MOLECULE_TYPE, label="Molekül Tipi")
    tree_type = forms.ChoiceField(choices=ALGORITMA, label="Filogenetik Ağaç Tipi")
    alignment_filetype = forms.ChoiceField(choices=ALIGNMENT_FILE_TYPE, label="Alignment Dosya Tipi")


class MultipleSequenceAlignmentModelForm(forms.ModelForm):
    class Meta:
        model = MultipleSequenceAlignment
        fields = ['in_file', 'method', 'alignment_filetype', 'molecule_type', 'tree_type']


class MaximumLikeHoodForm(forms.Form):
    file = forms.FileField(label="Fasta Dosyası Seçiniz")
    palm_tools = forms.ChoiceField(choices=PALM_TOOLS, label="PAML Aracı")
    molecule_type = forms.ChoiceField(choices=MOLECULE_TYPE, label="Molekül Tipi")
    tree_type = forms.ChoiceField(choices=ALGORITMA, label="Filogenetik Ağaç Tipi")






