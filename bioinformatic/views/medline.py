from django.shortcuts import *

from Bio import Entrez
from bioinformatic.models import MedlineArticle
from bioinformatic.forms.medline import MedlineArticleForm


def medline_article(request):
    form = MedlineArticleForm(request.POST or None)
    if request.method == "POST":
        if form.is_valid():
            form.save(commit=False)
            email = form.cleaned_data["email"]
            term = form.cleaned_data['term']

            Entrez.email = f"{email}"

            handle = Entrez.esearch(db="pubmed", term=term)

            record = Entrez.read(handle)

            idlist = record["IdList"]

            if len(idlist) == 0:
                return render(request, "bioinformatic/fasta/notfound.html",
                              {"msg": "Aradığınız Terim bulunamadı", "url": reverse("bioinformatic:medline_article"),
                               "bre": "Bulunamadı"})

            handle = Entrez.efetch(db="pubmed", id=idlist, rettype="medline", retmode="xml")

            records = Entrez.read(handle)

            if MedlineArticle.objects.exists():
                MedlineArticle.objects.all().delete()
                for record in records["PubmedArticle"]:
                    MedlineArticle.objects.create(title=record["MedlineCitation"]["Article"]["ArticleTitle"],
                                                  article_id=record["MedlineCitation"]["PMID"])
            else:
                for record in records["PubmedArticle"]:
                    MedlineArticle.objects.create(title=record["MedlineCitation"]["Article"]["ArticleTitle"],
                                                  article_id=record["MedlineCitation"]["PMID"])

            return render(request, "bioinformatic/medline/result.html",
                          {"articles": MedlineArticle.objects.all(), "bre": "Makaleler",
                           "count": MedlineArticle.objects.all().count()})

    return render(request, "bioinformatic/medline/search.html", {"form": form, "bre": "Makale Arama"})
