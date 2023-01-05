from django.contrib import admin
from .models import Profile, SocialMedia, ContactModel


@admin.register(ContactModel)
class ContactModelAdmin(admin.ModelAdmin):
    list_display = ['author', 'sender', 'contact_email', 'created']
    list_filter = ['author', 'created']
    search_fields = ['author']


@admin.register(Profile)
class ProfileAdmin(admin.ModelAdmin):
    list_display = ['user', 'birth_day', 'created']
    list_filter = ['user']
    search_fields = ['user']
    raw_id_fields = ['user']


@admin.register(SocialMedia)
class SocialMediaAdmin(admin.ModelAdmin):
    list_display = ['user', 'social', 'url']
    list_filter = ['user', 'social', 'url']
    search_fields = ['user']
    raw_id_fields = ['user']
