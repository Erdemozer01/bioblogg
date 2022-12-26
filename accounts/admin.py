from django.contrib import admin
from .models import Profile, Job, Education, SocialMedia


@admin.register(Profile)
class ProfileAdmin(admin.ModelAdmin):
    list_display = ['user', 'birth_day', 'created']
    list_filter = ['user']
    search_fields = ['user']
    raw_id_fields = ['user']


@admin.register(Job)
class JobAdmin(admin.ModelAdmin):
    list_display = ['user', 'job', 'department']
    list_filter = ['user', 'job', 'department']
    search_fields = ['user', 'job', 'department']
    raw_id_fields = ['user']


@admin.register(Education)
class EducationAdmin(admin.ModelAdmin):
    list_display = ['user', 'school', 'school_name', 'start', 'finish']
    list_filter = ['user', 'school', 'school_name']
    search_fields = ['user', 'school', 'school_name']


@admin.register(SocialMedia)
class SocialMediaAdmin(admin.ModelAdmin):
    list_display = ['user', 'social', 'url']
    list_filter = ['user', 'social', 'url']
    search_fields = ['user']
    raw_id_fields = ['user']
