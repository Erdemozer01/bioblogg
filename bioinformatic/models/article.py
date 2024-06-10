from django.db import models


class ArticleSearchModel(models.Model):
    user = models.ForeignKey('auth.User', on_delete=models.CASCADE, verbose_name='Kullanıcı', blank=True)

    email = models.EmailField()

    title = models.CharField(max_length=1000)

    article_id = models.CharField(max_length=100)

    created = models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi')

    def __str__(self):
        return self.title

    class Meta:
        db_table = 'articles'
        verbose_name = "Makale"
        verbose_name_plural = "Makale"
