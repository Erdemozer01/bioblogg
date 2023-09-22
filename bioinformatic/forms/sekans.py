from django import forms


class DNASekansForm(forms.Form):
    dna = forms.CharField(label="", widget=forms.Textarea(attrs={'placeholder': 'DNA SEKANSI GİRİNİZ'}))


class SequenceSlicingForm(forms.Form):
    sequence = forms.CharField(label="SEKANS", widget=forms.Textarea(attrs={'placeholder': 'SEKANS GİRİNİZ'}))
    kmer_len = forms.IntegerField(label="Kmer uzunluğu", min_value=1)
    dot = forms.IntegerField(label="Nükleotit pozisyonu", min_value=1)
    max_nuc = forms.IntegerField(label="Max Nükleotit pozisyonu", help_text="Nükleotit pozisyonundaki nükletit sayısı")
