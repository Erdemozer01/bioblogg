from django import forms

FILE_FORMAT = (
    ('', '-' * 30),
    ('fasta', 'FASTA'),
    ('genbank', 'GENBANK'),
)

PLOT_TYPE = (
    ('', '-' * 30),
    ('gc', "%GC Plot"),
    ('histogram', "Histogram Plot"),
    ('dot', "Dot Plot"),
    ('volcano', "Volcano Plot"),
)


class PlotSelectForm(forms.Form):
    plot_type = forms.ChoiceField(label="Plot Türü", choices=PLOT_TYPE)


class PlotForm(forms.Form):
    file_format = forms.ChoiceField(label="Dosya formatı seçiniz", choices=FILE_FORMAT)
    file = forms.FileField(label="Dosya seçiniz")


class VolcanoPlotForm(forms.Form):
    file = forms.FileField(label="CSV Dosya seçiniz",help_text=".* csv uzantılı olmalıdır")
