from django.db import models
from django.contrib.auth.models import User


class PasswordModel(models.Model):
    name = models.CharField(max_length=100, unique=True)
    password = models.CharField(max_length=2000, unique=True)
    created = models.DateTimeField(auto_now=True, verbose_name='Oluşturulma Tarihi', blank=True)

    def __str__(self):
        return self.name

    class Meta:
        verbose_name = 'Password'
        verbose_name_plural = 'Password'


def avatar(instance, filename):
    return 'users/avatar/{0}/{1}'.format(instance.user.profile, filename)


def cover(instance, filename):
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
    avatar = models.ImageField(upload_to=avatar, blank=True, verbose_name='Avatar:')
    gender = models.CharField(max_length=30, choices=Gender, default="", verbose_name='Cinsiyet', blank=True)
    location = models.CharField(max_length=100, verbose_name='Yaşadığı şehir', blank=True)
    about = models.TextField(verbose_name='Hakkımda', blank=True)
    job = models.CharField(max_length=100, verbose_name='Meslek', blank=True)
    phone = models.CharField(max_length=130, verbose_name='Telefon', blank=True, unique=True, null=True)
    birth_day = models.DateField(blank=True, null=True, verbose_name='Doğum Tarihi:')
    skills = models.CharField(max_length=2000, editable=True, verbose_name='Yetenekler', blank=True)
    web_site = models.URLField(verbose_name='Web Sitesi', blank=True)
    language = models.CharField(max_length=2000, editable=True, verbose_name='Bildiğiniz Diller', blank=True,
                                help_text="ingilizce")

    created = models.DateTimeField(auto_now=True, verbose_name='Katılma Tarihi', blank=True)

    def __str__(self):
        return str(self.user)

    class Meta:
        db_table = 'profile'
        verbose_name = 'Profil'
        verbose_name_plural = 'Profil'


class SocialMedia(models.Model):
    user = models.ForeignKey(Profile, on_delete=models.CASCADE, verbose_name='Kullanıcı Adı',
                             related_name='user_social')
    social = models.CharField(verbose_name="Sosyal Medya Adı", max_length=100, help_text="Facebook")
    url = models.URLField(verbose_name="Sosyal Medya Url Adresi", max_length=100)
    created = models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi', blank=True)

    def __str__(self):
        return self.social

    class Meta:
        db_table = 'profile_social'
        verbose_name = 'Sosyal Medya Bilgileri'
        verbose_name_plural = 'Sosyal Medya Bilgileri'


class ContactModel(models.Model):
    title = models.CharField(max_length=150, verbose_name="Konu:", blank=True)
    receiver = models.ForeignKey(User, on_delete=models.CASCADE, verbose_name='Alıcı: ',
                                 related_name='messages_receiver')
    sender = models.ForeignKey('auth.User', on_delete=models.CASCADE, verbose_name='Gönderen', related_name='gönderen')
    content = models.TextField(verbose_name='Mesaj')
    contact_email = models.EmailField()
    is_read = models.BooleanField(default=False, verbose_name='Okundu')
    is_report = models.BooleanField(default=False, verbose_name='Kötüye kullanım')
    created = models.DateTimeField(auto_now_add=True, verbose_name='Gönderilme Tarihi', blank=True)

    def __str__(self):
        return self.sender.username

    class Meta:
        db_table = 'contact'
        verbose_name = 'Mesajlar'
        verbose_name_plural = 'Mesajlar'
        ordering = ['-created']
