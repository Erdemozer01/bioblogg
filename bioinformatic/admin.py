from django.contrib import admin
from .models import FileModel, BioinformaticModel, RecordModel

admin.AdminSite.login_template = "registration/login.html"
admin.AdminSite.site_header = "BioBlog"


class CreateDatabaseModelInline(admin.StackedInline):
    model = RecordModel
    extra = 1


class FileModelModelInline(admin.StackedInline):
    model = FileModel
    extra = 1


@admin.register(BioinformaticModel)
class BioinformaticModelAdmin(admin.ModelAdmin):
    list_display = ['user', 'tool', "molecule", 'reading_file_format', 'writing_file_format', 'created']
    list_filter = ['tool', 'created']
    search_fields = ['tool', 'reading_file_format', 'writing_file_format']
    inlines = [
        FileModelModelInline,
        CreateDatabaseModelInline
    ]
