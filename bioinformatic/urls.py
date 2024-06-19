from django.urls import path
from bioinformatic import views

app_name = "bioinformatic"

urlpatterns = [
    # Sekans
    path('anasayfa/', views.bioinformatic_home, name="home"),
    path('dna-sekans-okuması/', views.SekansView.sequence_analiz, name="dna_seq_read"),
    path('dna-sekans-kesme/', views.SekansView.Kmer_SeqSlicing, name="dna_seq_slice"),
    path('dna-sekans-translate/', views.SekansView.translation, name="dna_seq_translate"),
    path('create_frame_seq/', views.SekansView.create_frame_seq, name="create_frame_seq"),
    path('alignment-score/', views.SekansView.alignment_score, name="alignment_score"),
    path('temperature-melting/', views.SekansView.TemperatureMeltingView, name="temp_melt"),

    # Reading
    path('file-reading/<user>/', views.FileReadingView.file_reading, name="file_reading"),
    path('phylo-tree-creating/', views.FileReadingView.PhylogeneticTree, name="pyhlo_tree"),
    path('blast/', views.FileReadingView.blast, name="blast"),
    path('3d-molekul-goruntuleme/', views.FileReadingView.molecule_viewer, name="3d_molecule_view"),

    # Writing
    path('file-format-select/<user>/', views.FileWritingView.file_writing_format_select,
         name="file_writing_format_select"),
    path('dosya-oluşturma/<format>/<user>/', views.FileWritingView.FileWritingView, name="file_writing"),
    path('file-writing-list/<format>/<user>/', views.FileWritingListView.as_view(), name="file_writing_list"),
    path('download/<user>/<format>/', views.CreateFileView, name='create_and_download'),
    path('record-detail/<pk>/<user>/<description>/', views.RecordDetailView.as_view(), name='record_detail'),
    path('delete/<pk>/', views.RecordDeleteView, name='delete_record'),

    # Download
    path('download/', views.download_file, name='download'),

    # FileConvert
    path('dosya-dönüştürme/', views.FileConvert, name='file_convert'),

    # Entrez Tools
    path('entrez-tools/', views.entrez_tools, name='entrez_tools'),

    # Entrez Tools
    path('molecule-demo/', views.molecule_demo, name='molecule_demo'),

]
