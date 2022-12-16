from django.contrib import admin
from .models import Category, Posts, Comments


@admin.register(Comments)
class CommentsAdmin(admin.ModelAdmin):
    list_display = ['commentator', 'post', 'created']


@admin.register(Category)
class CategoryAdmin(admin.ModelAdmin):
    list_display = ['title', 'created']


@admin.register(Posts)
class PostAdmin(admin.ModelAdmin):
    list_display = ['title', 'author', 'category', 'created']
    list_filter = ['author', 'category', 'created']
    search_fields = ['title', 'author']
    search_help_text = "Başlık ve yazara göre arama"
