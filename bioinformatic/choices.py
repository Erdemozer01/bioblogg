BLAST_PROGRAM = (
    ('', '------------'),
    ('blastn', 'BLASTN'),
    ('blastp', 'BLASTP'),
    ('blastx', 'BLASTX'),
    ('tblastn', 'TBLASTN'),
    ('tblastn', 'TBLASTN'),
)

BLAST_DATABASE = (
    ('', '------------'),
    ('nr', 'nr'),
    ('nt', 'nt'),
)

ENTREZ_SELECT = (
    ('', '------------'),
    ('art', 'GÜNCEL MAKALE ARAMA'),
    ('nuc', 'NUCLEOTIDE'),
)

RECORD_FORMAT = (
    ('', '------------'),
    ('file', 'DOSYA'),
    ('gi', 'Gİ(GENİNFO) Numarası'),
)

TOOLS = (
    ('', '------------'),
    ('SEKANS İŞLEMLERİ', 'SEKANS İŞLEMLERİ'),
    ('DOSYA İŞLEMLERİ', 'DOSYA İŞLEMLERİ'),
)

TRANS_TABLE = (
    ("", "------------"),
    ("1", "Standart Kod"),
    ("2", "Omurgalı Mitokondri Kodu"),
    ("3", "Maya Mitokondri Kodu"),
    ("4", "Küf, Protozoon ve Kölenterat Mitokondri Kodu ve Mikoplazma/Spiroplasma Kodu"),
    ("5", "Omurgasız Mitokondri Kodu"),
    ("6", "Siliat, Dasikladas ve Heksamita Nükleer Kodu"),
    ("9", "Ekinoderm ve Yassı Solucan Mitokondri Kodu"),
    ("10", "Euplotid Nükleer Kodu"),
    ("11", "Bakteriyel, Arkeal ve Bitki Plastid Kodu, prokaryotik virüsler"),
    ("12", "Alternatif Maya Nükleer Kodu"),
    ("13", "Ascidian Mitokondri Kodu"),
    ("14", "Alternatif Yassı Solucan Mitokondri Kodu"),
    ("16", "Klorofis Mitokondri Kodu"),
    ("21", "Trematod Mitokondriyal Kodu"),
    ("22", "Scenedesmus obliquus Mitokondri Kodu"),
    ("23", "Thraustochytrium Mitokondri Kodu"),
    ("24", "Rhabdopleuridae Mitokondri Kodu"),
    ("25", "Aday Bölüm SR1 ve Gracilibacteria Kodu"),
    ("26", "Pachysolen tannophilus Nükleer Kodu"),
    ("27", "Karyorelict Nükleer Kodu"),
    ("28", "Kondilostoma Nükleer Kodu"),
    ("29", "Mezodinyum Nükleer Kodu"),
    ("30", "Peritrich Nükleer Kodu"),
    ("31", "Blastocrithidia Nükleer Kodu"),
    ("33", "Cephalodiscidae Mitokondriyal UAA-Tyr Kodu"),
)

READ_FILE_FORMAT = [
    ("abi", "ABİ"),
    ("abi-trim", "abi-trim".upper()),
    ("ace", "ace".upper()),
    ("cif-atom", "cif-atom".upper()),
    ("cif-seqres", "cif-seqres".upper()),
    ("clustal", "clustal".upper()),
    ("embl", "embl".upper()),
    ("embl-cds", "embl-cds".upper()),
    ("fasta", "fasta".upper()),
    ("fastq", "fastq".upper()),
    ("fasta-2line", "fasta-2line".upper()),
    ("fastq-sanger", "fastq-sanger".upper()),
    ("fastq-solexa", "fastq-solexa".upper()),
    ("fastq-illumina", "fastq-illumina".upper()),
    ("genbank", "GENBANK"),
    ("gck", "gck".upper()),
    ("ig", "ig".upper()),
    ("imgt", "imgt".upper()),
    ("nib", "nib".upper()),
    ("nexus", "nexus".upper()),
    ("pdb-atom", "pdb-atom".upper()),
    ("phd", "phd".upper()),
    ("phylip", "phylip".upper()),
    ("pir", "pir".upper()),
    ("qual", "qual".upper()),
    ("seqxml", "seqxml".upper()),
    ("sff", "sff".upper()),
    ("sff-trim", "sff-trim".upper()),
    ("snapgene", "snapgene".upper()),
    ("stockholm", "stockholm".upper()),
    ("swiss", "swiss".upper()),
    ("tab", "tab".upper()),
    ("qual", "qual".upper()),
    ("twobit", "twobit".upper()),
    ("uniprot-xml", "uniprot-xml".upper()),
    ("xdna", "xdna".upper()),
]

WRITE_FILE_FORMAT = [
    ("clustal", "clustal".upper()),
    ("embl", "embl".upper()),
    ("fasta", "fasta".upper()),
    ("fastq", "fastq".upper()),
    ("fasta-2line", "fasta-2line".upper()),
    ("fastq-sanger", "fastq-sanger".upper()),
    ("fastq-solexa", "fastq-solexa".upper()),
    ("fastq-illumina", "fastq-illumina".upper()),
    ("gb", "GENBANK (.gb) "),
    ("genbank", "GENBANK (.genbank)"),
    ("imgt", "imgt".upper()),
    ("nexus", "nexus".upper()),
    ("phd", "phd".upper()),
    ("phylip", "phylip".upper()),
    ("pir", "pir".upper()),
    ("seqxml", "seqxml".upper()),
    ("sff", "sff".upper()),
    ("stockholm", "stockholm".upper()),
    ("tab", "tab".upper()),
    ("qual", "qual".upper()),
    ("xdna", "xdna".upper()),
]

ALIGNMENT_MODE = (
    ('', '------------'),
    ('local', 'LOCAL'),
    ('global', 'GLOBAL'),
)

MATRIS = (
    ('', '------------'),
    ('BENNER22', 'BENNER22'),
    ('BENNER6', 'BENNER6'),
    ('BENNER74', 'BENNER74'),
    ('BLOSUM45', 'BLOSUM45'),
    ('BLOSUM50', 'BLOSUM50'),
    ('BLOSUM62', 'BLOSUM62'),
    ('BLOSUM80', 'BLOSUM80'),
    ('BLOSUM90', 'BLOSUM90'),
    ('DAYHOFF', 'DAYHOFF'),
    ('FENG', 'FENG'),
    ('HOXD70', 'HOXD70'),
    ('JOHNSON', 'JOHNSON'),
    ('JONES', 'JONES'),
    ('LEVIN', 'LEVIN'),
    ('MCLACHLAN', 'MCLACHLAN'),
    ('MDM78', 'MDM78'),
    ('PAM250', 'PAM250'),
    ('PAM30', 'PAM30'),
    ('PAM70', 'PAM70'),
    ('RAO', 'RAO'),
    ('RISLER', 'RISLER'),
    ('SCHNEIDER', 'SCHNEIDER'),
    ('STR', 'STR'),
    ('TRANS', 'TRANS'),
)

METHOD = (
    ('', '------------'),
    ('muscle', 'MUSCLE'),
    ('clustalw2', 'CLUSTALW2'),
    ('omega', 'ClustalOmega'),
    ('paml', 'Maximum Likelihood (PAML)'),
)

TREE_ALGORITMA = (
    ('', '------------'),
    ('nj', 'Neighbor Joining'),
    ('upgma', 'UPGMA'),
)

MOLECULE_TYPE = (
    ('', '------------'),
    ('DNA', 'DNA'),
    ('RNA', 'RNA'),
    ('protein', 'PROTEİN'),
)

ALIGNMENT_FILE_TYPE = (
    ('', '------------'),
    ('fasta', 'FASTA'),
    ('clustal', 'CLUSTAL'),
    ('phylip', 'PHYLİB'),
    ('nexus', 'NEXUS'),
)

PALM_TOOLS = (
    ('', '------------'),
    ('baseml', 'BASEML'),
    ('basemlg', 'BASEMLG'),
    ('codeml', 'CODEML'),
)
