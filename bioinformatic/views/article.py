from Bio import Entrez
from bioinformatic.forms import ArticleForm
from django.shortcuts import *
from django.contrib import messages
from bioinformatic.models.article import ArticleSearchModel


def pubmed(request):
    if request.user.is_anonymous:
        from django.conf import settings
        messages.error(request, "Lütfen Giriş Yapınız")
        return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))

    form = ArticleForm(request.POST or None)

    if request.method == "POST":
        if ArticleSearchModel.objects.filter(user=request.user).exists():
            ArticleSearchModel.objects.filter(user=request.user).delete()

        if form.is_valid():
            email = form.cleaned_data["email"]
            term = form.cleaned_data['term']

            Entrez.email = email

            handle = Entrez.esearch(db="pubmed", term=term)

            record = Entrez.read(handle)

            idlist = record["IdList"]

            if len(idlist) == 0:
                return render(request, "exception/page-404.html",
                              {"msg": "Aradığınız Terim bulunamadı", "url": reverse("bioinformatic:article"),
                               "bre": "Bulunamadı"})

            handle = Entrez.efetch(db="pubmed", id=idlist, rettype="medline", retmode="xml")

            records = Entrez.read(handle)

            for record in records["PubmedArticle"]:
                ArticleSearchModel.objects.create(
                    user=request.user,
                    title=record["MedlineCitation"]["Article"]["ArticleTitle"],
                    article_id=record["MedlineCitation"]["PMID"],
                )

            articles = ArticleSearchModel.objects.filter(user=request.user)

            return render(request, "bioinformatic/articles.html", {"articles": articles})

    return render(request, "bioinformatic/form.html", {'form': form, 'title': "Güncel Makale Arama"})
