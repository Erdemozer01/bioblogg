# Generated by Django 4.1.4 on 2023-01-05 11:23

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('accounts', '0002_alter_contactmodel_email'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='job',
            name='user',
        ),
        migrations.DeleteModel(
            name='Education',
        ),
        migrations.DeleteModel(
            name='Job',
        ),
    ]
