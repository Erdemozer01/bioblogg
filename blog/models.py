import time
from django.db import models
from django.utils import timezone
from autoslug.fields import AutoSlugField
from hitcount.models import HitCount, Hit
from django.contrib.contenttypes.fields import GenericRelation
from django_ckeditor_5.fields import CKEditor5Field, CKEditor5Widget

num_time = time.time()


def cover_upload_to(instance, filename):
    local_time = time.localtime(num_time)
    return 'blog/posts/{username}/{year}-{mon}-{day}/{hour}_{minute}/{filename}'.format(
        username=instance.author.username,
        year=local_time.tm_year,
        mon=local_time.tm_mon,
        day=local_time.tm_mday,
        hour=local_time.tm_hour,
        minute=local_time.tm_min,
        filename=filename
    )


def extra_img_upload_to(instance, filename):
    local_time = time.localtime(num_time)
    return 'blog/posts/{username}/{year}-{mon}-{day}/{hour}_{minute}/extra_image/{filename}'.format(
        username=instance.author.username,
        year=local_time.tm_year,
        mon=local_time.tm_mon,
        day=local_time.tm_mday,
        hour=local_time.tm_hour,
        minute=local_time.tm_min,
        filename=filename
    )


class Category(models.Model):
    image = models.ImageField(upload_to='blog/category/', verbose_name='Kategori Fotosu:')
    title = models.CharField(max_length=100, verbose_name='Kategori:')
    explain = CKEditor5Field(verbose_name='Kategori Tanımı:', config_name='extends')
    slug = AutoSlugField(populate_from='title', unique=True)

    publish = models.DateTimeField(default=timezone.now, verbose_name='Yayınlama Tarihi')
    created = models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi')
    updated = models.DateTimeField(auto_now=True, verbose_name='Güncellenme Tarihi')

    def __str__(self):
        return self.title

    class Meta:
        db_table = 'category'
        verbose_name = 'Kategori'
        verbose_name_plural = 'Kategori'


class PublishedManager(models.Manager):
    def get_queryset(self):
        return super().get_queryset() \
            .filter(status=Posts.Status.PUBLISHED)


class Posts(models.Model):
    class Status(models.TextChoices):
        DRAFT = 'DF', 'Taslak'
        PUBLISHED = 'PB', 'Yayınla'

    category = models.ForeignKey(Category, verbose_name="Kategori", related_name='post', on_delete=models.CASCADE)
    cover = models.ImageField(upload_to=cover_upload_to, verbose_name="Gönderi Fotosu:")
    title = models.CharField(max_length=250, verbose_name="Başlık:", blank=False)
    text = CKEditor5Field(verbose_name='İçerik', config_name='extends')
    slug = AutoSlugField(populate_from="title", unique=True)
    author = models.ForeignKey('auth.User', on_delete=models.CASCADE, verbose_name='Yazar')

    tags = models.CharField(max_length=255, blank=True, verbose_name="Etiketler")
    hit_count = GenericRelation(HitCount, object_id_field='object_pk', related_query_name='hit_count_generic_relation')
    likes = models.ManyToManyField('auth.User', related_name='post_like', verbose_name="Beğendim", blank=True)
    dislike = models.ManyToManyField('auth.User', related_name='post_dislike', verbose_name="Beğenmedim", blank=True)

    extra_image1 = models.ImageField(upload_to=extra_img_upload_to, verbose_name="Slayt1",
                                     blank=True)
    extra_image2 = models.ImageField(upload_to=extra_img_upload_to, verbose_name="Slayt2",
                                     blank=True)
    extra_image3 = models.ImageField(upload_to=extra_img_upload_to, verbose_name="Slayt3",
                                     blank=True)

    created = models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi')
    updated = models.DateTimeField(auto_now=True, verbose_name='Güncellenme Tarihi')

    status = models.CharField(max_length=2, choices=Status.choices, default=Status.DRAFT,
                              verbose_name='Yayınlanma Durumu')
    objects = models.Manager()
    published = PublishedManager()

    def __str__(self):
        return self.title

    class Meta:
        db_table = 'posts'
        ordering = ['created']
        indexes = [
            models.Index(fields=['created']),
        ]
        verbose_name_plural = 'Gönderi'
        verbose_name = 'Gönderi'


class ReadManager(models.Manager):
    def get_queryset(self):
        return super().get_queryset() \
            .filter(status=Comments.Status.UNREAD)


class Comments(models.Model):
    class Status(models.TextChoices):
        READ = 'okundu', 'Okundu'
        UNREAD = 'okunmadı', 'okunmadı'

    commentator = models.ForeignKey('auth.User', on_delete=models.CASCADE, verbose_name='Yazan:', related_name="post")
    post = models.ForeignKey(Posts, on_delete=models.CASCADE, verbose_name='Gönderi:', related_name="comment")
    comment = models.TextField(verbose_name='Yorum Yap:')

    report = models.ManyToManyField('auth.User', related_name='comment_report', verbose_name="Kötüye Kullanım",
                                    blank=True)
    likes = models.ManyToManyField('auth.User', related_name='comment_likes', verbose_name="Beğendim", blank=True)
    dislike = models.ManyToManyField('auth.User', related_name='comment_dislike', verbose_name="Beğenmedim", blank=True)
    status = models.CharField(max_length=30, choices=Status.choices, default=Status.UNREAD,
                              verbose_name='Okunma Durumu')

    publish = models.DateTimeField(default=timezone.now, verbose_name='Yayınlama Tarihi')
    created = models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi')
    updated = models.DateTimeField(auto_now=True, verbose_name='Güncellenme Tarihi')

    objects = models.Manager()  # The default manager.
    published = ReadManager()  # Our custom manager.

    def __str__(self):
        return str(self.commentator) + ',' + self.post.title

    class Meta:
        db_table = 'comments'
        verbose_name = 'Yorum'
        verbose_name_plural = 'Yorumlar'


class Subscribe(models.Model):
    email = models.EmailField()

    def __str__(self):
        return self.email

    class Meta:
        db_table = 'subscribe'
        verbose_name = 'abone'
        verbose_name_plural = 'aboneler'


class BlogContactModel(models.Model):
    map_url = models.TextField(verbose_name="Google map url")
    title = models.CharField(max_length=250, verbose_name="Başlık")
    explain = models.TextField(verbose_name="Açıklama")
    address = models.CharField(max_length=100, verbose_name="Adres")
    telephone = models.CharField(max_length=100, verbose_name="Telefon")
    email = models.EmailField(verbose_name="Email adresimiz")
    created = models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi')

    def __str__(self):
        return self.title

    class Meta:
        db_table = 'blog_contact'
        verbose_name = 'İletişim'
        verbose_name_plural = 'İletişim'


class BlogContactFormModel(models.Model):
    name = models.CharField(max_length=100, verbose_name="Ad Soyad")
    email = models.CharField(verbose_name="Email", max_length=100)
    message = models.TextField(verbose_name="Mesaj")
    created = models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi')

    def __str__(self):
        return self.name

    class Meta:
        db_table = "blog_contact_form_model"
        verbose_name = "Blog İletişim"
        verbose_name_plural = "Blog İletişim"


class About(models.Model):
    image = models.ImageField(upload_to='blog/posts/about/')
    title = models.CharField(max_length=100, verbose_name="Başlık")
    text = models.TextField(verbose_name="Hakkımızda")
