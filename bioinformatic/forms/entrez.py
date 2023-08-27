from django import forms

DATABASE = (
    ("pubmed", "Pubmed"),
    ("protein", "Protein"),
    ("nucleotide", "Nucleotide"),
    ("genome", "Genome"),
    ("gene", "Gene"),
    ("structure", "Structure"),
    ("taxonomy", "Taxonomy"),
    ("ipg", "İpg"),
    ("nuccore", "Nuccore"),
    ("annotinfo", "Annotinfo"),
    ("assembly", "Assembly"),
    ("bioproject", "Bioproject"),
    ("biosample", "Biosample"),
    ("blastdbinfo", "Blastdbinfo"),
    ("books", "Books"),
    ("cdd", "cdd"),
    ("clinvar", "Clinvar"),
    ("gap", "GAP"),
    ("gapplus", "Gapplus"),
    ("grasp", "Grasp"),
    ("grasp", "Grasp"),
    ("dbvar", "dbvar"),
    ("gds", "GDS (GEO datasets)"),
    ("geoprofiles", "Geoprofiles"),
    ("homologene", "Homologene"),
    ("medgen", "Medgen"),
    ("mesh", "mesh"),
    ("nlmcatalog", "nlmcatalog"),
    ("omim", "omim"),
    ("orgtrack", "orgtrack"),
    ("pmc", "pmc"),
    ("popset", "popset"),
    ("proteinclusters", "Proteinclusters"),
    ("pcassay", "pcassay"),
    ("protfam", "protfam"),
    ("pccompound", "pccompound"),
    ("seqannot", "seqannot"),
    ("biocollections", "biocollections"),
    ("snp", "snp"),
    ("sra", "sra"),
    ("gtr", "gtr"),
)

RETTYPE = (
    ("gb", "Genbank"),
    ("fasta", "Fasta"),
    ("xml", "XML"),
    ("medline", "Medline"),
)


class EntrezForm(forms.Form):
    email = forms.EmailField()

    database = forms.ChoiceField(choices=DATABASE, label="Veritabanı Seçiniz",
                                 widget=forms.Select(attrs={'placeholder': "Veri Tabanı Seçiniz"}))

    accession = forms.CharField(label="Erişim yada id Numarası")

    rettype = forms.ChoiceField(label="Dosya Türü", choices=RETTYPE)
