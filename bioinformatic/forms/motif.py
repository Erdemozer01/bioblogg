from django import forms


class JasparMotifCreateForm(forms.Form):
    file = forms.FileField(label="DNA Fasta Dosyası Giriniz")
    seq_length = forms.IntegerField(label="Motif Sekans Uzunluğu")
    data_type = forms.CharField(label="Veri Tipi")
    name = forms.CharField(label="Jaspar Adı")
    matrix_id = forms.CharField(label="Matrix İD")
    tax_group = forms.CharField(label="Taxonomik Grup")

class FastaDNASeqMotifForm(forms.Form):
    file = forms.FileField(label="DNA Fasta Dosyası Giriniz")
    length = forms.IntegerField(label="Motif Sekans Uzunluğu")
