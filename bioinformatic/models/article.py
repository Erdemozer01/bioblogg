from django.db import models


class ArticleSearchModel(models.Model):
    user = models.ForeignKey('auth.User', on_delete=models.CASCADE, verbose_name='Kullanıcı', blank=True)

    title = models.CharField(max_length=1000, verbose_name='Başlık')

    article_id = models.CharField(max_length=100, verbose_name="Makale ID")

    created = models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi')

    def __str__(self):
        return self.title

    class Meta:
        db_table = 'articles'
        verbose_name = "Güncel Makale Arama"
        verbose_name_plural = "Güncel Makale Arama"
