from django.db import models
from pathlib import Path
from bioinformatic.choices import *

BASE_DIR = Path(__file__).resolve().parent.parent


def read_file(instance, filename):
    return 'laboratory/reading/{username}/{username}_{filename}'.format(
        username=instance.user.username, filename=filename)


def alignment_file(instance, filename):
    return 'laboratory/alignment/{username}/{username}_{filename}'.format(
        username=instance.user.username, filename=filename)


class BioinformaticAnalizModel(models.Model):
    user = models.ForeignKey('auth.User', on_delete=models.CASCADE, verbose_name='Laborant', blank=True)
    tool = models.CharField(max_length=100, verbose_name="İşlem", blank=True)
    reading_file_format = models.CharField(choices=READ_FILE_FORMAT, max_length=100, verbose_name="Dosya Okuma Formatı",
                                           blank=True)
    writing_file_format = models.CharField(choices=WRITE_FILE_FORMAT, max_length=100,
                                           verbose_name="Dosya Yazdırma Formatı", blank=True)
    molecule = models.CharField(choices=MOLECULE_TYPE, max_length=100, verbose_name="Molekül", blank=True)
    read_file = models.FileField(verbose_name="Dosya", upload_to=read_file, blank=True,
                                 help_text="Not : Max. dosya boyutu 5mb olmalıdır.")
    trans_table = models.CharField(choices=TRANS_TABLE, max_length=100, verbose_name="Dönüşüm Tablosu", blank=True)
    to_stop = models.BooleanField(verbose_name="Stop kodonları dahil et.", default=False)

    molecule_id = models.CharField(verbose_name="İD", max_length=100, blank=True)
    name = models.CharField(verbose_name="Adı", max_length=100, blank=True)
    description = models.CharField(verbose_name="TANIMLAMA", max_length=100, blank=True)
    dbxrefs = models.CharField(verbose_name="dbxrefs", max_length=100, blank=True)
    source = models.CharField(verbose_name="source", max_length=100, blank=True)
    accession = models.CharField(verbose_name="accession", max_length=100, blank=True)
    author = models.CharField(max_length=100, verbose_name="YAZAR", blank=True, null=True)
    keywords = models.CharField(max_length=100, verbose_name="Anahtar Kelimeler", blank=True, null=True)
    organism = models.CharField(verbose_name="ORGANİZMA", max_length=100, blank=True, null=True)
    taxonomy = models.CharField(verbose_name="TAKSONOMİ", max_length=100, blank=True)
    strain = models.CharField(max_length=100, blank=True)
    gene = models.CharField(max_length=100, blank=True)
    select = models.BooleanField(default=False)
    email = models.EmailField(blank=True)
    seq = models.TextField(verbose_name="SEKANS")
    seq_len = models.IntegerField(verbose_name="Sekans Uzunluğu", blank=True, null=True)
    gc = models.FloatField(verbose_name="%GC", blank=True, null=True)
    pro_seq = models.TextField(verbose_name="PROTEİN SEKANS", blank=True)
    pro_seq_len = models.IntegerField(verbose_name="PROTEİN SEKANS UZUNLUĞU", blank=True, null=True)
    alignment = models.TextField(verbose_name="Alignment", blank=True)
    alignment_file = models.FileField(verbose_name="Alignment Dosyası", upload_to=alignment_file, blank=True)
    out_file_fasta = models.FileField(verbose_name="Sonuç Fasta", upload_to="laboratory/file/fasta/", blank=True)
    out_file_genbank = models.FileField(verbose_name="Sonuç Genbank", upload_to="laboratory/file/genbank/", blank=True)
    out_file_xml = models.FileField(verbose_name="Sonuç XML", upload_to="laboratory/file/xml/", blank=True)
    out_file_tree = models.FileField(verbose_name="Sonuç TREE", upload_to="laboratory/file/tree/", blank=True)
    created = models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi')

    def __str__(self):
        return str(self.user)

    class Meta:
        db_table = "bioinformatic"
        verbose_name = "Biyoinformatik Analiz"
        verbose_name_plural = "Biyoinformatik Analiz"


class BioinformaticDatabaseModel(models.Model):
    user = models.ForeignKey('auth.User', on_delete=models.CASCADE, verbose_name='Ekleyen (Kullanıcı ADI)', blank=True)
    db_name = models.CharField(max_length=100, blank=True, verbose_name="Veri Tabanı Adı")
    format = models.CharField(choices=READ_FILE_FORMAT, max_length=100, verbose_name="Format", blank=True)
    molecule = models.CharField(choices=MOLECULE_TYPE, max_length=100, verbose_name="Molekül", blank=True)
    created = models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi')

    def __str__(self):
        return str(self.user) + "-" + self.db_name

    class Meta:
        verbose_name = "Biyoinformatik Veri Tabanı"
        verbose_name_plural = "Biyoinformatik Veri Tabanı"


class CreateDatabaseModel(models.Model):
    record = models.ForeignKey(BioinformaticDatabaseModel, on_delete=models.CASCADE, blank=True)
    name = models.CharField(max_length=100, blank=True)
    description = models.CharField(max_length=100, blank=True)
    organism = models.CharField(max_length=100, blank=True)
    taxonomy = models.CharField(max_length=100, blank=True)
    annotations = models.CharField(max_length=100, blank=True)
    letter_annotations = models.CharField(max_length=100, blank=True)
    dbxrefs = models.CharField(max_length=100, blank=True)
    features = models.CharField(max_length=100, blank=True)
    author = models.CharField(max_length=100, blank=True)
    keywords = models.CharField(max_length=100, blank=True)
    sequence = models.TextField(blank=True)

    class Meta:
        verbose_name = "Biyoinformatik Kayıt"
        verbose_name_plural = "Biyoinformatik Kayıt"
