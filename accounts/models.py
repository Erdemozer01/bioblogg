from django.db import models
from django.contrib.auth.models import User


def avatar(instance, filename):
    return 'users/avatar/{0}/{1}'.format(instance.user.profile, filename)


def cover(instance, filename):
    # file will be uploaded to MEDIA_ROOT/beat/author/<filename>
    return 'users/cover/{0}/{1}'.format(instance.user.profile, filename)


class Profile(models.Model):
    Gender = [
        ('', '------'),
        ('Erkek', 'Erkek'),
        ('Kadın', 'Kadın'),
        ('Belirtmek istemiyorum', 'Belirtmek istemiyorum'),
    ]
    user = models.OneToOneField(User, on_delete=models.CASCADE, verbose_name='Kullanıcı Adı', related_name='profile')
    first_name = models.CharField(max_length=150, verbose_name="Ad")
    last_name = models.CharField(max_length=150, verbose_name="Soyad")
    email = models.EmailField(verbose_name="Email")
    cover = models.ImageField(upload_to=cover, blank=True, verbose_name='Kapak Fotosu:')
    avatar = models.ImageField(upload_to=avatar, blank=True, verbose_name='Avatar:')
    gender = models.CharField(max_length=30, choices=Gender, default="", verbose_name='Cinsiyet', blank=True)
    location = models.CharField(max_length=100, verbose_name='Yaşadığı şehir', blank=True)
    about = models.TextField(verbose_name='Hakkımda', blank=True)
    job = models.CharField(max_length=100, verbose_name='Meslek', blank=True)
    phone = models.CharField(max_length=13, verbose_name='Telefon', blank=True, unique=True, null=True)
    birth_day = models.DateField(blank=True, null=True, verbose_name='Doğum Tarihi:')
    skils = models.CharField(max_length=200, editable=True, verbose_name='Yetenekler', blank=True)

    facebook = models.URLField(verbose_name='Facebook', blank=True)
    twitter = models.URLField(verbose_name='Twitter', blank=True)
    instagram = models.URLField(verbose_name='İnstagram', blank=True)

    created = models.DateTimeField(auto_now=True, verbose_name='Katılma Tarihi', blank=True)

    def __str__(self):
        return str(self.user)

    class Meta:
        db_table = 'profile'
        verbose_name = 'Profil'
        verbose_name_plural = 'Profil'
