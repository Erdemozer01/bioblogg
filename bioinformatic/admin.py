from django.contrib import admin
from .models import BioinformaticAnalizModel, BioinformaticDatabaseModel, CreateDatabaseModel

admin.AdminSite.login_template = "registration/login.html"
admin.AdminSite.site_header = "BioBlog"


@admin.register(BioinformaticAnalizModel)
class BioinformaticModelAdmin(admin.ModelAdmin):
    list_display = ['user', 'tool', "molecule", 'reading_file_format', 'created']
    list_filter = ['tool', 'created']
    search_fields = ['tool', 'file_format']


class BioinformaticDatabaseModelInline(admin.StackedInline):
    model = CreateDatabaseModel
    extra = 1


@admin.register(BioinformaticDatabaseModel)
class CreateDatabaseModelAdmin(admin.ModelAdmin):
    list_display = ["user", "db_name", "format", "molecule", "created"]
    list_filter = ["db_name", "format", "molecule", "created"]
    search_fields = ["db_name", "format"]
    inlines = [
        BioinformaticDatabaseModelInline
    ]
