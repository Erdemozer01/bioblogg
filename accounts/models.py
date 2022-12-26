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
    language = models.CharField(max_length=200, editable=True, verbose_name='Bildiğiniz Diller', blank=True, help_text="ingilizce")

    created = models.DateTimeField(auto_now=True, verbose_name='Katılma Tarihi', blank=True)

    def __str__(self):
        return str(self.user)

    class Meta:
        db_table = 'profile'
        verbose_name = 'Profil'
        verbose_name_plural = 'Profil'


class Education(models.Model):
    user = models.ForeignKey(User, on_delete=models.CASCADE, verbose_name='Kullanıcı Adı', related_name='user_edu')
    school = models.CharField(verbose_name="Okul Türü", max_length=100, help_text="Üniversite")
    school_name = models.CharField(verbose_name="Okul", max_length=100, help_text="Akdeniz Üniversitesi")
    start = models.CharField(verbose_name="Başlangıç yılı", max_length=100, help_text="2010")
    finish = models.CharField(verbose_name="Bitiş yılı:", max_length=100, help_text="2014")
    explain = models.TextField(verbose_name='Açıklama:', blank=True)
    created = models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi', blank=True)

    def __str__(self):
        return self.school + "|" + self.school_name

    class Meta:
        db_table = 'education'
        verbose_name = 'Eğitim Bilgileri'
        verbose_name_plural = 'Eğitim Bilgileri'


class Job(models.Model):
    user = models.ForeignKey(User, on_delete=models.CASCADE, verbose_name='Kullanıcı Adı', related_name='user_job')
    job = models.CharField(verbose_name="Meslek", max_length=100, help_text="Mühendis")
    department = models.CharField(verbose_name="Pozisyon", max_length=100, help_text="Gıda Mühendisi")
    location = models.CharField(verbose_name="İş yeri adı", max_length=100, help_text="Sabancı Holding")
    start = models.CharField(verbose_name="Başlangıç yılı", max_length=100, help_text="2010")
    explain = models.TextField(verbose_name='Açıklama:', blank=True)
    created = models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi', blank=True)

    def __str__(self):
        return self.job + "|" + self.department

    class Meta:
        db_table = 'job'
        verbose_name = 'İş Bilgileri'
        verbose_name_plural = 'İş Bilgileri'


class SocialMedia(models.Model):
    user = models.ForeignKey(User, on_delete=models.CASCADE, verbose_name='Kullanıcı Adı', related_name='user_social')
    social = models.CharField(verbose_name="Sosyal Medya Adı", max_length=100, help_text="Facebook", unique=True)
    url = models.URLField(verbose_name="Sosyal Medya Url Adresi", max_length=100)
    created = models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi', blank=True)

    def __str__(self):
        return self.social

    class Meta:
        db_table = 'profile_social'
        verbose_name = 'Sosyal Medya Bilgileri'
        verbose_name_plural = 'Sosyal Medya Bilgileri'
