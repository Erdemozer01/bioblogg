from django.db import models
from pathlib import Path
from bioinformatic.choices import *

BASE_DIR = Path(__file__).resolve().parent.parent


def file_upload_to(instance, filename):
    return 'laboratory/{username}/{filename}'.format(
        username=instance.records.user, filename=filename)


class BioinformaticModel(models.Model):
    user = models.ForeignKey('auth.User', on_delete=models.CASCADE, verbose_name='İşlemi yapan', blank=True)
    tool = models.CharField(max_length=100, verbose_name="İşlem", blank=True)
    reading_file_format = models.CharField(choices=READ_FILE_FORMAT, max_length=100, verbose_name="Dosya Okuma Formatı",
                                           blank=True)
    writing_file_format = models.CharField(choices=WRITE_FILE_FORMAT, max_length=100,
                                           verbose_name="Dosya Yazdırma Formatı", blank=True)
    molecule = models.CharField(choices=MOLECULE_TYPE, max_length=100, verbose_name="Molekül", blank=True)
    trans_table = models.CharField(choices=TRANS_TABLE, max_length=100, verbose_name="Dönüşüm Tablosu", blank=True)
    to_stop = models.BooleanField(verbose_name="Stop kodonları dahil et.", default=False)
    db_name = models.CharField(max_length=100, blank=True, verbose_name="Veri Tabanı Adı")
    created = models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi')

    def __str__(self):
        return str(self.user) + "-" + self.db_name

    class Meta:
        verbose_name = "Biyoinformatik"
        verbose_name_plural = "Biyoinformatik"


class FileModel(models.Model):
    records = models.ForeignKey(BioinformaticModel, on_delete=models.CASCADE, blank=True, related_name="records_files")
    file = models.FileField(verbose_name="Dosya", upload_to=file_upload_to, blank=True,
                            help_text="Not : Max. dosya boyutu 5mb olmalıdır.")

    class Meta:
        verbose_name = "Dosya"
        verbose_name_plural = "Dosyalar"


class RecordModel(models.Model):
    records = models.ForeignKey(BioinformaticModel, on_delete=models.CASCADE, blank=True, related_name="record_content")

    molecule = models.CharField(max_length=100, blank=True, verbose_name="Molekül")
    molecule_id = models.CharField(max_length=100, blank=True, verbose_name="Molekül İD")
    name = models.CharField(max_length=100, blank=True, verbose_name="Molekül adı")
    description = models.CharField(max_length=100, blank=True, verbose_name="Tanımlama")
    letter_annotations = models.CharField(max_length=100, blank=True, null=True)
    db_xrefs = models.CharField(max_length=100, blank=True, null=True, verbose_name="veritabanı referans")
    features = models.TextField(blank=True, null=True)
    gene = models.CharField(max_length=100, blank=True, null=True)

    annotations = models.CharField(max_length=100, blank=True)
    taxonomy = models.CharField(max_length=100, blank=True)
    date = models.CharField(max_length=100, blank=True)
    data_file_division = models.CharField(max_length=100, blank=True)
    source = models.CharField(max_length=100, blank=True, null=True)
    topology = models.CharField(max_length=100, choices=TOPOLOGY, blank=True)
    accession = models.CharField(max_length=100, blank=True, null=True)
    keywords = models.CharField(max_length=100, blank=True, null=True)
    organism = models.CharField(max_length=100, blank=True, null=True)

    sequence = models.TextField(verbose_name="Sekans", blank=True)
    seq_len = models.IntegerField(verbose_name="Sekans Uzunluğu", blank=True, null=True)
    protein_sequence = models.TextField(verbose_name="Protein Sekans", blank=True)
    pro_seq_len = models.IntegerField(verbose_name="PROTEİN SEKANS UZUNLUĞU", blank=True, null=True)
    gc = models.FloatField(verbose_name="%GC", blank=True, null=True)

    class Meta:
        verbose_name = "Biyoinformatik Kayıt"
        verbose_name_plural = "Biyoinformatik Kayıtlar"
