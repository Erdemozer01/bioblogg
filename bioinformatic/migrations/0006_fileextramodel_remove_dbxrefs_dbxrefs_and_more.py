# Generated by Django 4.1.6 on 2023-09-17 12:19

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('bioinformatic', '0005_alter_letter_annotations_letter_annotations'),
    ]

    operations = [
        migrations.CreateModel(
            name='FileExtraModel',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('dbxrefs', models.CharField(blank=True, max_length=100, verbose_name='dbxrefs')),
                ('features', models.CharField(blank=True, max_length=100, verbose_name='features')),
                ('annotations', models.CharField(blank=True, max_length=100, verbose_name='annotations adı')),
                ('annotations_description', models.CharField(max_length=100, verbose_name='annotations')),
                ('letter_annotations', models.CharField(blank=True, max_length=100, verbose_name='letter_annotations adı')),
                ('letter_annotations_description', models.CharField(max_length=100, verbose_name='letter_annotations')),
            ],
        ),
        migrations.RemoveField(
            model_name='dbxrefs',
            name='dbxRefs',
        ),
        migrations.RemoveField(
            model_name='features',
            name='features',
        ),
        migrations.RemoveField(
            model_name='letter_annotations',
            name='letter_annotations',
        ),
        migrations.AlterField(
            model_name='bioinformaticanalizmodel',
            name='seq',
            field=models.TextField(verbose_name='SEKANS'),
        ),
        migrations.DeleteModel(
            name='Annotations',
        ),
        migrations.DeleteModel(
            name='DBXRefs',
        ),
        migrations.DeleteModel(
            name='Features',
        ),
        migrations.DeleteModel(
            name='Letter_Annotations',
        ),
        migrations.AddField(
            model_name='fileextramodel',
            name='parent_model',
            field=models.ForeignKey(blank=True, on_delete=django.db.models.deletion.CASCADE, related_name='file_extra_information', to='bioinformatic.bioinformaticanalizmodel'),
        ),
    ]