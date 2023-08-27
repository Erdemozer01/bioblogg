import os.path
from django.core.files.uploadedfile import InMemoryUploadedFile
from django.contrib.auth.models import User
from django.db import models
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent


class LabSlideModel(models.Model):
    image = models.ImageField(upload_to="laboratory/slide/")
    title = models.CharField(max_length=100, verbose_name='Başlık')
    text = models.CharField(max_length=100, verbose_name='Alt Başlık:')
    content = models.TextField(verbose_name='İçerik')

    def __str__(self):
        return self.title + '|' + self.subtitle

    class Meta:
        db_table = "bioinformatic_home"
        verbose_name = "Slide"
        verbose_name_plural = "Slide"


class FastaRead(models.Model):
    gene = models.CharField(max_length=1000)
    sequence = models.TextField()

    def __str__(self):
        return self.gene

    class Meta:
        db_table = "fasta_read"
        verbose_name = "Fasta Okuması"
        verbose_name_plural = "Fasta Okumaları"


class GenbankRead(models.Model):
    organism = models.CharField(verbose_name="Organizma", max_length=1000)
    description = models.CharField(max_length=1000, verbose_name="Tanım")
    protein_id = models.CharField(max_length=255, verbose_name="protein_id", blank=True, null=True)
    taxonomy = models.CharField(max_length=1000, verbose_name="Taksonomi", null=True, blank=True)
    dna_sequence = models.TextField()
    dna_sequence_len = models.BigIntegerField(verbose_name="DNA Uzunluğu", null=True, blank=True)
    protein_sequence = models.TextField(blank=True, null=True)
    protein_sequence_len = models.BigIntegerField(verbose_name="Protein Uzunluğu", null=True, blank=True)

    def __str__(self):
        return self.organism

    class Meta:
        ordering = ["id"]
        db_table = "genbank_read"
        verbose_name = "Genbank Okuması"
        verbose_name_plural = "Genbank Okumaları"


class PubMedArticle(models.Model):
    email = models.EmailField()

    article_id = models.CharField(max_length=100)

    title = models.CharField(max_length=1000)

    link = models.URLField(default="https://pubmed.ncbi.nlm.nih.gov/")

    created = models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi')

    def __str__(self):
        return self.link

    class Meta:
        verbose_name = "PubMed Makale"
        verbose_name_plural = "PubMed Makale"


class MedlineArticle(models.Model):
    email = models.EmailField()

    article_id = models.CharField(max_length=1000)

    doi = models.CharField(max_length=1000, default="")

    title = models.CharField(max_length=1000)

    term = models.CharField(max_length=1000)

    created = models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi')

    def __str__(self):
        return self.title


class SwissProtModel(models.Model):
    accessions = models.CharField(max_length=1000, verbose_name="Erişim Numarası")
    organism = models.CharField(max_length=1000, verbose_name="Organizma")
    sequence = models.TextField(verbose_name="Sekans")
    sequence_length = models.CharField(max_length=1000, verbose_name="Sekans Uzunluğu")

    def __str__(self):
        return self.accessions

    class Meta:
        db_table = "swiss-prot"
        verbose_name = "swiss-prot"
        verbose_name_plural = "swiss-prot"


class BigFileUploadModel(models.Model):
    big_file = models.FileField(verbose_name="Dosya Seçme", upload_to="laboratory/bigfile/")
    created = models.DateTimeField(auto_now_add=True, null=True)

    def __str__(self):
        return self.big_file

    class Meta:
        db_table = "bigfile"
        verbose_name = "bigfile"
        verbose_name_plural = "bigfile"


class FileUploadModel(models.Model):
    file = models.FileField(verbose_name="Dosya Seçme", upload_to="laboratory/bigfile/")
    created = models.DateTimeField(auto_now_add=True, null=True)

    def __str__(self):
        return self.file

    class Meta:
        db_table = "file"
        verbose_name = "file"
        verbose_name_plural = "files"


MSA_TYPE = (
    ('', '------------'),
    ('muscle', 'Muscle'),
    ('clustalw2', 'ClustalW2'),
    ('omega', 'Clustal Omega'),
    ('paml', 'Maximum Likelihood (PAML)'),
)

TREE_TYPE = (
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


def upload_to(instance, filename):
    return 'MultipleSequenceAlignment/{username}/{username}_{filename}'.format(
        username=instance.user.username, filename=filename)


def molecule_upload_to(instance, filename):
    return 'molecule/{username}/{filename}'.format(
        username=instance.user.username, filename=filename.lower())


def molecule_upload_to2(instance, filename):
    return 'molecule/{username}/{filename}'.format(
        username=instance.user.username, filename=filename)


def molecule_path():
    return os.path.join(BASE_DIR, "media", "molecule")


class MolecularModel(models.Model):
    user = models.ForeignKey('auth.User', on_delete=models.CASCADE, verbose_name='Laborant')
    if InMemoryUploadedFile:
        in_file = models.FileField(verbose_name="PDB Dosyası", upload_to=molecule_upload_to2, blank=True, null=True)
    else:
        in_file = models.FileField(verbose_name="PDB Dosyası", upload_to=molecule_upload_to, blank=True, null=True)
    file_path = models.FilePathField(path=molecule_path, blank=True, null=True, recursive=True)
    file_name = models.CharField(max_length=100, verbose_name="Dosya Adı:", blank=True, null=True)
    name = models.CharField(max_length=100, verbose_name="Adı:", blank=True, null=True)
    author = models.CharField(max_length=100, verbose_name="Yazar:", blank=True, null=True)
    head = models.CharField(max_length=100, verbose_name="Özellik:", blank=True, null=True)
    id_code = models.CharField(max_length=100, verbose_name="id cod:", blank=True, null=True)
    keywords = models.CharField(max_length=100, verbose_name="Kategori:", blank=True, null=True)
    id_name = models.CharField(max_length=100, verbose_name="id:", blank=True, null=True)
    created = models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi')

    def __str__(self):
        return str(self.user)


class MultipleSequenceAlignment(models.Model):
    user = models.ForeignKey('auth.User', on_delete=models.CASCADE, verbose_name='Laborant')
    method = models.CharField(choices=MSA_TYPE, verbose_name="Multiple Sekans Alignment Aracı", max_length=1000)
    tree_type = models.CharField(choices=TREE_TYPE, verbose_name="Filogenetik Ağaç Türü", max_length=1000)
    molecule_type = models.CharField(choices=MOLECULE_TYPE, verbose_name="Molekül Tipi", max_length=1000)
    palm_tools = models.CharField(choices=PALM_TOOLS, verbose_name="Maximum Likelihood (PAML)", max_length=1000,
                                  blank=True, null=True)
    alignment_filetype = models.CharField(choices=ALIGNMENT_FILE_TYPE, verbose_name="Alignment Dosya Tipi",
                                          max_length=1000)
    in_file = models.FileField(verbose_name="Girdi Dosyası", upload_to=upload_to, blank=True, null=True)
    out_file = models.FileField(verbose_name="Çıktı Dosyası", upload_to=upload_to, blank=True, null=True)
    align_file = models.FileField(verbose_name="Aligned Dosyası", upload_to=upload_to, blank=True, null=True)
    stats = models.FileField(verbose_name="İstatistik Dosyası", upload_to=upload_to, blank=True, null=True)
    scores = models.FileField(verbose_name="Alignment Scor Dosyası", upload_to=upload_to, blank=True, null=True)
    ml_file = models.FileField(verbose_name="Maximum Likelihood", upload_to=upload_to, blank=True, null=True)
    tree_file = models.FileField(verbose_name="Filogenetik Ağaç Dosyası", upload_to=upload_to, blank=True, null=True)
    alignment_chart = models.FileField(verbose_name="Aligment Haritası", upload_to=upload_to, blank=True, null=True)
    tree = models.ImageField(verbose_name="Filogenetik Ağaç", blank=True, null=True)
    cluster_csv = models.FileField(verbose_name="ClusterGram CSV Dosyası", upload_to=upload_to, blank=True, null=True)
    created = models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi')

    def __str__(self):
        return self.method

    class Meta:
        db_table = "msa"
        verbose_name = "Multiple Sequence Alignment"
        verbose_name_plural = "Multiple Sequence Alignment"


class BiologicalResourcesDatabases(models.Model):
    name = models.CharField(max_length=100, verbose_name="Kaynak")
    url = models.URLField(verbose_name="URL")
    created = models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi')

    def __str__(self):
        return self.name

    class Meta:
        verbose_name = "Biyolojik Kaynaklar ve Veri Tabanları"
        verbose_name_plural = "Biyolojik Kaynaklar ve Veri Tabanları"


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
    ('nr', 'BLASTN'),
    ('nt', 'BLASTP'),
)


def upload_to_blast(instance, filename):
    return 'blast/{username}/{username}_{filename}'.format(
        username=instance.user.username, filename=filename)


class BlastQueryResults(models.Model):
    user = models.ForeignKey('auth.User', on_delete=models.CASCADE, verbose_name='Laborant')
    program = models.CharField(choices=BLAST_PROGRAM, verbose_name="Blast Programı Seçiniz", max_length=1000)
    database = models.CharField(choices=BLAST_DATABASE, verbose_name="Blast Programı Seçiniz", max_length=1000)
    input_file = models.FileField(verbose_name="Fasta Dosyası", upload_to=upload_to_blast)
    output_file = models.FileField(verbose_name="BlastXML Dosyası", upload_to=upload_to_blast, blank=True, null=True)
    hsp_file = models.FileField(verbose_name="HSP Dosyası", upload_to=upload_to_blast, blank=True, null=True)
    created = models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi')

    def __str__(self):
        return self.user

    class Meta:
        verbose_name = "Blast Metodu"
        verbose_name_plural = "Blast Metodu"


def upload_to_motif(instance, filename):
    return 'motif/{username}/{username}_{filename}'.format(
        username=instance.user.username, filename=filename)


class FastaDNAMotifModel(models.Model):
    user = models.ForeignKey('auth.User', on_delete=models.CASCADE, verbose_name='Laborant')
    motif_file = models.FileField(verbose_name="Motif Nucleotit Sekans Dosyası", upload_to=upload_to_motif)
    mpf = models.FileField(verbose_name="Nucleotide matrix position", upload_to=upload_to_motif, blank=True,
                           null=True)
    pwm = models.FileField(verbose_name="Compute position weight matrices", upload_to=upload_to_motif, blank=True,
                           null=True)

    pssm = models.FileField(verbose_name="Compute position weight matrices", upload_to=upload_to_motif, blank=True,
                            null=True)
    created = models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi')

    def __str__(self):
        return str(self.user)

    class Meta:
        verbose_name = "Fasta DNA Sekans Motif"
        verbose_name_plural = "Fasta DNA Sekans Motif"


from reportlab.lib import colors

COLORS = (
    ('', '------------'),
    ((colors.blue), 'Mavi'),
    ("black", 'SİYAH'),
    ("red", 'KIRMIZI'),
    ("green", 'YEŞİL'),
    ("gold", 'SARI'),
    ("orange", 'TURUNCU'),
)


def upload_to_diagram(instance, filename):
    return 'diagram/{username}/{username}_{filename}'.format(
        username=instance.user.username, filename=filename)


class FileFormat(models.Model):
    name = models.CharField(max_length=255, verbose_name="Dosya Formatı")

    def __str__(self):
        return self.name


DIAGRAM_SHAPE = (
    ('', '-' * 30),
    ('ARROW', 'OK'),
    ('BOX', 'KUTU'),
)

DIAGRAM_FORMAT = (
    ('', '------------'),
    ('linear', 'DÜZ'),
    ('circular', 'DAİRESEL'),
)


class DiagramModel(models.Model):
    user = models.ForeignKey('auth.User', on_delete=models.CASCADE, verbose_name='Laborant')
    out_file = models.FileField(blank=True, null=True, upload_to=upload_to_diagram)
    image = models.ImageField(upload_to=upload_to_diagram)

    def __str__(self):
        return str(self.user)


class RestrictionModel(models.Model):
    user = models.ForeignKey(DiagramModel, on_delete=models.CASCADE, verbose_name='Laborant')
    enzymes = models.CharField(max_length=1000, verbose_name="Enzim Adı:", blank=True, null=True)
    site = models.CharField(max_length=1000, verbose_name="Bölge:", blank=True, null=True)

    def __str__(self):
        return self.enzymes


def upload_to_graphic(instance, filename):
    return 'graphic/{username}/{username}_{filename}'.format(
        username=instance.user.username, filename=filename)


class GraphicModels(models.Model):
    user = models.ForeignKey('auth.User', on_delete=models.CASCADE, verbose_name='Laborant')
    graph_type = models.CharField(max_length=1000, verbose_name="Grafik Türü")
    format = models.CharField(max_length=1000, verbose_name="Dosya Formatı", blank=True, null=True)
    histogram_plot = models.ImageField(upload_to=upload_to_graphic, verbose_name="Histogram", blank=True, null=True)
    gc_plot = models.ImageField(upload_to=upload_to_graphic, verbose_name="%GC Plot", blank=True, null=True)
    dot_plot = models.ImageField(upload_to=upload_to_graphic, verbose_name="Dot Plot", blank=True, null=True)
    created = models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi')

    def __str__(self):
        return str(self.user)


class MaximumFileSize(models.Model):
    file_size = models.BigIntegerField()
    created = models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi')

    def __str__(self):
        return str(self.file_size)
