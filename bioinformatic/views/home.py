from django.views import generic
from bioinformatic.models import LabSlideModel


class BioinformaticHomeView(generic.ListView):
    template_name = "bioinformatic/home.html"
    model = LabSlideModel

    def get_context_data(self, *, object_list=None, **kwargs):
        context = super(BioinformaticHomeView, self).get_context_data(**kwargs)
        context['bre'] = "Biyoinformatik Anasayfa"
        return context
