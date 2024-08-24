from django.urls import path
from bioinformatic import views

app_name = "bioinformatic"

urlpatterns = [
    # Sekans
    path('anasayfa/', views.bioinformatic_home, name="home"),
    path('dna-sekans-okuması/', views.SekansView.sequence_analiz, name="dna_seq_read"),
    path('kmer-oluşturma/', views.SekansView.Kmer_SeqSlicing, name="dna_seq_slice"),
    path('dna-sekans-translate/', views.SekansView.translation, name="dna_seq_translate"),
    path('frame-seq/', views.create_frame_seq, name="create_frame_seq"),
    path('alignment-score/', views.SekansView.alignment_score, name="alignment_score"),
    path('temperature-melting/', views.SekansView.TemperatureMeltingView, name="temp_melt"),

    # Reading
    path('file-reading/', views.FileReadingView.file_reading, name="file_reading"),

    path('blast/', views.FileReadingView.blast, name="blast"),
    path('alignment-mapping/', views.alignment_mapping, name="alignment_mapping"),


    # Writing
    path('file-format-select/<user>/', views.FileWritingView.file_writing_format_select,
         name="file_writing_format_select"),

    path('dosya-oluşturma/<user>/<format>/', views.FileWritingView.FileWritingView, name="file_writing"),

    path('olustur/<user>/<format>/', views.CreateFileView, name='create_file'),

    path('delete/<pk>/', views.RecordDeleteView, name='file_create_rec_delete'),

    # Download
    path('download/', views.download_file, name='download_file'),

    # FileConvert
    path('dosya-dönüştürme/', views.FileConvert, name='file_convert'),

    # Makale
    path('entrez/', views.EntrezView, name='entrez_view'),

    # molecule_viewer
    path('single-3d-molekul-goruntuleme/', views.single_molecule_view, name='single_molecule_3d_view'),
    path('multiple-3d-molekul-goruntuleme/', views.multi_molecule_view, name='multiple_molecule_3d_view'),
    path('molecule-2d-görüntüleme/', views.molecule_2d_view, name='molecule_2d_view'),




]
