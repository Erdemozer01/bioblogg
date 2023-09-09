from django.db import models


class BiologicalResourcesDatabases(models.Model):
    name = models.CharField(max_length=100, verbose_name="Kaynak")
    url = models.URLField(verbose_name="URL")
    created = models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi')

    def __str__(self):
        return self.name

    class Meta:
        verbose_name = "Biyolojik Kaynaklar ve Veri Tabanları"
        verbose_name_plural = "Biyolojik Kaynaklar ve Veri Tabanları"
