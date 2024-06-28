from .home import bioinformatic_home
from .SekansView import sequence_analiz, translation, alignment_score, create_frame_seq
from .FileReadingView import file_reading, PhylogeneticTree
from .FileWritingView import file_writing_format_select, CreateFileView, RecordDetailView, RecordDeleteView, \
    FileWritingListView
from .FileDownloadView import download_file
from .FileConvertView import FileConvert
from .Entrez import EntrezToolsView
from .molecule_viewer import single_molecule_view, multi_molecule_view
