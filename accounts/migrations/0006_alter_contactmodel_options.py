# Generated by Django 4.1.4 on 2023-01-05 21:51

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('accounts', '0005_alter_contactmodel_options_alter_contactmodel_table'),
    ]

    operations = [
        migrations.AlterModelOptions(
            name='contactmodel',
            options={'verbose_name': 'Mesajlar', 'verbose_name_plural': 'Mesajlar'},
        ),
    ]
