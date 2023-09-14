from django.urls import path
from bioinformatic import views

app_name = "bioinformatic"

urlpatterns = [

    path('anasayfa/', views.BioinformaticHomeView.as_view(), name="home"),
    path('dna-sekans-okuması/', views.SekansView.sequence_analiz, name="dna_seq_read"),
    path('dna-sekans-translate/', views.SekansView.translation, name="dna_seq_translate"),
    # Reading
    path('file-reading/<user>/', views.FileReadingView.file_reading, name="file_reading"),
    path('file-reading/<user>/sonuçlar/', views.FileReadingView.FileReadingResultView.as_view(),
         name="file_reading_results"),
    path('file-reading-detail/<slug:user>/<pk>/<slug:description>/', views.FileReadingView.FileReadDetailView.as_view(),
         name="file_reading_detail"),
    path('file-reading/<user>/protein/', views.FileReadingView.ProteinPickView, name="file_reading_pick_protein"),
    path('istatistik-verileri/<user>/', views.FileReadingView.stats_view, name="stats"),
    path('alignment-score/<user>/', views.FileReadingView.alignment_score, name="alignment_score"),
    # Writing

]
