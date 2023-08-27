from django import forms

BLAST_PROGRAM = (
    ('', '------------'),
    ('blastn', 'BLASTN'),
    ('blastp', 'BLASTP'),
    ('blastx', 'BLASTX'),
    ('tblastn', 'TBLASTN'),
    ('tblastn', 'TBLASTN'),
)

BLAST_DATABASE = (
    ('', '------------'),
    ('nr', 'nr'),
    ('nt', 'nt'),
)


class BlastXMLForm(forms.Form):
    input_file = forms.FileField(label="Fasta Dosyası")
    program = forms.ChoiceField(choices=BLAST_PROGRAM, label="Blast Programı Seçiniz")
    database = forms.ChoiceField(choices=BLAST_DATABASE, label="Blast Veritabanı Seçiniz")

