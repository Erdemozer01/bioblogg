import time
from django.db import models
from django.utils import timezone
from ckeditor_uploader.fields import RichTextUploadingField
from autoslug.fields import AutoSlugField
from hitcount.models import HitCount, Hit
from django.contrib.contenttypes.fields import GenericRelation

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
    explain = models.TextField(verbose_name='Kategori Tanımı:')
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

    cover = models.ImageField(upload_to=cover_upload_to, verbose_name="Gönderi Fotosu:")
    title = models.CharField(max_length=250, verbose_name="Başlık:", blank=False)
    text = RichTextUploadingField(verbose_name='İçerik', blank=False)
    slug = AutoSlugField(populate_from="title", unique=True)
    author = models.ForeignKey('auth.User', on_delete=models.CASCADE, verbose_name='Yazar')
    category = models.ForeignKey(Category, verbose_name="Kategori", related_name='post', on_delete=models.CASCADE)
    tags = models.CharField(max_length=255, blank=True, verbose_name="Etiketler")
    hit_count = GenericRelation(HitCount, object_id_field='object_pk', related_query_name='hit_count_generic_relation')

    extra_image1 = models.ImageField(upload_to=extra_img_upload_to, verbose_name="Gönderiye Foto Ekleme",
                                     blank=True)
    extra_image2 = models.ImageField(upload_to=extra_img_upload_to, verbose_name="Gönderiye Foto Ekleme",
                                     blank=True)
    extra_image3 = models.ImageField(upload_to=extra_img_upload_to, verbose_name="Gönderiye Foto Ekleme",
                                     blank=True)

    created = models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi')
    updated = models.DateTimeField(auto_now=True, verbose_name='Güncellenme Tarihi')

    status = models.CharField(max_length=2, choices=Status.choices, default=Status.DRAFT,
                              verbose_name='Yayınlanma Durumu')
    objects = models.Manager()  # The default manager.
    published = PublishedManager()  # Our custom manager.

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


class Comments(models.Model):
    commentator = models.ForeignKey('auth.User', on_delete=models.CASCADE, verbose_name='Yazan:', related_name="post")
    post = models.ForeignKey(Posts, on_delete=models.CASCADE, verbose_name='Gönderi:', related_name="comment")
    comment = models.TextField(verbose_name='Yorum Yap:')

    report = models.ManyToManyField('auth.User', related_name='comment_report', verbose_name="Kötüye Kullanım")
    likes = models.ManyToManyField('auth.User', related_name='comment_likes', verbose_name="Beğendim")
    dislike = models.ManyToManyField('auth.User', related_name='comment_dislike', verbose_name="Beğenmedim")

    publish = models.DateTimeField(default=timezone.now, verbose_name='Yayınlama Tarihi')
    created = models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi')
    updated = models.DateTimeField(auto_now=True, verbose_name='Güncellenme Tarihi')

    def __str__(self):
        return str(self.commentator) + ', ' + self.post.title[:40]

    class Meta:
        db_table = 'comments'
        verbose_name = 'Yorum'
        verbose_name_plural = 'Yorumlar'
