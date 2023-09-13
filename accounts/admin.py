from django.contrib import admin
from .models import Profile, SocialMedia, ContactModel, PasswordModel


@admin.register(PasswordModel)
class PasswordModelAdmin(admin.ModelAdmin):
    list_display = ['name', 'created']


@admin.register(ContactModel)
class ContactModelAdmin(admin.ModelAdmin):
    list_display = ['sender', 'receiver', 'contact_email', 'is_read', 'created']
    list_filter = ['sender', 'created']
    search_fields = ['sender']


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
