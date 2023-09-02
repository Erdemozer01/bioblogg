from django.db import models


class Notifications(models.Model):
    user = models.ForeignKey("auth.User", on_delete=models.CASCADE, verbose_name='Kullanıcı Adı')
    type = models.CharField(max_length=100, verbose_name="Bildirim türü")
    title = models.CharField(max_length=100, verbose_name="Bildirim başlığı")
    text = models.TextField(verbose_name="Bildirim içeriği", blank=True)
    url = models.URLField(blank=True)
    is_read = models.BooleanField(default=False)
    created = models.DateTimeField(auto_now=True, verbose_name='Oluşturulma Tarihi')

    def __str__(self):
        return str(self.user)

    class Meta:
        db_table = "notifications"
        verbose_name = "Bildrim"
        verbose_name_plural = "Bildrimler"
