from django.urls import path
from bioinformatic import views

app_name = "bioinformatic"

urlpatterns = [
    # Sekans
    path('anasayfa/', views.BioinformaticHomeView.as_view(), name="home"),
    path('dna-sekans-okuması/', views.SekansView.sequence_analiz, name="dna_seq_read"),
    path('dna-sekans-kesme/', views.SekansView.SequenceSlicing, name="dna_seq_slice"),
    path('dna-sekans-translate/', views.SekansView.translation, name="dna_seq_translate"),

    # Reading
    path('file-reading/<user>/', views.FileReadingView.file_reading, name="file_reading"),
    path('file-reading/<user>/sonuçlar/', views.FileReadingView.FileReadingResultView.as_view(),
         name="file_reading_results"),
    path('file-reading-detail/<slug:user>/<pk>/', views.FileReadingView.FileReadDetailView.as_view(),
         name="file_reading_detail"),
    path('file-reading/<user>/protein/', views.FileReadingView.ProteinPickView, name="file_reading_pick_protein"),
    path('istatistik-verileri/<user>/', views.FileReadingView.stats_view, name="stats"),
    path('alignment-score/<user>/', views.FileReadingView.alignment_score, name="alignment_score"),

    # Writing
    path('file-format-select/<user>/', views.FileWritingView.file_writing_format_select,
         name="file_writing_format_select"),
    path('dosya-oluşturma/<format>/<user>/', views.FileWritingView.FileWritingView, name="file_writing"),
    path('file-writing-list/<format>/<user>/', views.FileWritingListView.as_view(), name="file_writing_list"),
    path('download/<user>/<format>/', views.CreateFileView, name='create_and_download'),
    path('record-detail/<pk>/<user>/<format>/', views.RecordDetailView.as_view(), name='record_detail'),
    path('delete/<pk>/', views.RecordDeleteView, name='delete_record'),

    # Download
    path('download/', views.download_file, name='download'),

    # FileConvert
    path('dosya-dönüştürme/', views.FileConvert, name='file_convert'),
]
