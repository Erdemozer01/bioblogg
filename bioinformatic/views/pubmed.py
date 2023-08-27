from django.shortcuts import *
from bioinformatic.forms.pubmed import PubMedForm
from Bio import Entrez
from bioinformatic.models import PubMedArticle

import requests
from bs4 import BeautifulSoup
import re


def pubmed(request):
    form = PubMedForm(request.POST or None)
    if request.method == "POST":
        if form.is_valid():

            article_id = form.cleaned_data['article_id']

            email = form.cleaned_data['email']

            Entrez.email = f"{email}"

            record = Entrez.read(Entrez.elink(dbfrom="pubmed", id=article_id))

            if PubMedArticle.objects.exists():
                PubMedArticle.objects.all().delete()
                for link in record[0]["LinkSetDb"][0]["Link"]:
                    PubMedArticle.objects.create(link="https://pubmed.ncbi.nlm.nih.gov/{}".format(link["Id"]),
                                                 article_id=link["Id"])
            else:
                for link in record[0]["LinkSetDb"][0]["Link"]:
                    PubMedArticle.objects.create(link="https://pubmed.ncbi.nlm.nih.gov/{}".format(link["Id"]),
                                                 article_id=link["Id"])

        return render(request, "bioinformatic/pubmed/result.html",
                      {'link': PubMedArticle.objects.all(), "bre": "Sonuçlar",
                       "count": PubMedArticle.objects.all().count()})

    return render(request, "bioinformatic/pubmed/search.html", {"form": form, "bre": "PubMed Id Makale Arama"})
