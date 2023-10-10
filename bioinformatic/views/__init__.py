from .home import BioinformaticHomeView
from .SekansView import sequence_analiz, translation, alignment_score
from .FileReadingView import (
    file_reading, FileReadingResultView, ProteinPickView, stats_view, PhylogeneticTree
)
from .FileWritingView import file_writing_format_select, CreateFileView, RecordDetailView, RecordDeleteView, \
    FileWritingListView
from .download import download_file
from .FileConvertView import FileConvert