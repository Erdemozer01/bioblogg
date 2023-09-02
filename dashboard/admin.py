from django.contrib import admin
from .models import Notifications


# Register your models here.

@admin.register(Notifications)
class NotificationsModelAdmin(admin.ModelAdmin):
    list_display = ['user', 'type', 'is_read', 'created']
    list_filter = ['user', 'type', 'created']
    search_fields = ['user', 'type']
