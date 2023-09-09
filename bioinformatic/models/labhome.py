from django.db import models
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent


class LabSlideModel(models.Model):
    image = models.ImageField(upload_to="laboratory/slide/")
    title = models.CharField(max_length=100, verbose_name='Başlık')
    text = models.CharField(max_length=100, verbose_name='Alt Başlık:')
    content = models.TextField(verbose_name='İçerik')

    def __str__(self):
        return self.title

    class Meta:
        db_table = "bioinformatic_home"
        verbose_name = "Slide"
        verbose_name_plural = "Slide"