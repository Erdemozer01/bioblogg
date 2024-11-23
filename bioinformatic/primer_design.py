from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils import MeltingTemp as mt

fasta = "ATGACAGAGGAGGAAATAACAGCTCTTGTGATCGATAATGGTTCTGGGATGGTGAAAGCGGGTTTCGCGGGAGACGACGCGCCCCGTGCAGTTTTCCCATCCATAGTGGGACGTCCGCGACACCAAGCAGTCATGGTTGGCATGGGGCAAAAAGACTCTTATGTGGGCGACGAAGCGCAGTCGAAACGCGGTATCCTATCGCTCAAGTACCCGATCGAGCACGGCATCGTAACGAACTGGGATGATATGGAAAAGATCTGGTACCACACTTTTTACAACGAGTTGCGCATCTCCCCCGAGGATCATCCAGTGTTGTTAACCGAAGCACCTTTGAACCCGAAAGCCAACAGGGAAAAAATGACCCAGATCATGTTCGAGACGTTCAACGTTCCCGCGATGTACGTCGCGATCCAGGCAGTATTGTCGTTGTATGCTTCGGGGCGCACCACGGGTATTGTCGTCGACAGCGGCGACGGTGTGACGCACACTGTGCCCATCTACGAAGGCTATGCCTTGCCGCATGCTATCATGCGCATAGATCTGGCTGGACGAGACCTAACCGACTATCTCGCCAAGATCCTTACAGAGCGCGGATATTCATTCACGACGACAGCCGAACGAGAGATCGTGCGAGACATCAAGGAAAAATGTTGTTATGTGGCACAGGACTACGATCATGAGTTGGAGATTGCCTCATCGCAGCCGGCGAAAATCGATAAGCAGTATGAACTCCCCGACGGACAAATCATCACGATTGGGAGTGAGCGTTTCAGATGCCCGGAGGTGCTCTTCCAGCCGTCCTTGATCGGTATGGAGGGCGAGGGTATTCATAATGTTGCTTATCAGAGCATCATGAAATGCGATGTCGACATCCGGAAAGATCTGTACGCAAACGTGGTTCTCAGCGGCGGCACGACGATGTTCCCAGGCATCGCGGATCGGATGCAGCGGGAACTCGCTAGTGTTGCACCATCTTCGGTGAAGATCAAGCTTGTAGCGCCAGCGGAGCGCAAATATAGCGTGTGGATCGGCGGCAGCATTTTGGCCTCGTTGAGCACTTTTCAGCAGATGTGGATTAGTAAGGCGGAGTATGACGAGTTCGGACCCTCTGTGGTACACCGCAAATGTTTCTGA"

def design_primers(dna_sequence, primer_length=20, gc_min=40, gc_max=60, tm_min=50, tm_max=60):
    """
    DNA dizisi için primerler tasarlar.

    Args:
        dna_sequence: DNA dizisi (str).
        primer_length: İstenen primer uzunluğu (int).
        gc_min: Minimum GC içeriği (%).
        gc_max: Maksimum GC içeriği (%).
        tm_min: Minimum erime sıcaklığı (°C).
        tm_max: Maksimum erime sıcaklığı (°C).

    Returns:
        Bir sözlük:
            - "forward_primer": İleri primer dizisi (str).
            - "reverse_primer": Ters primer dizisi (str).
            - "forward_gc": İleri primerin GC içeriği (float).
            - "reverse_gc": Ters primerin GC içeriği (float).
            - "forward_tm": İleri primerin erime sıcaklığı (float).
            - "reverse_tm": Ters primerin erime sıcaklığı (float).
        Uygun primer bulunamadığında None döner.
    """

    seq = Seq(dna_sequence)
    for i in range(len(seq) - 2 * primer_length + 1):
        print(i)
        forward = seq[i:i + primer_length]
        reverse = seq[i + primer_length:i + 2 * primer_length].reverse_complement()

        gc_f = gc_fraction(forward) * 100
        gc_r = gc_fraction(reverse) * 100
        tm_f = mt.Tm_NN(forward)  # Daha doğru sonuçlar için Tm_NN kullanıyoruz
        tm_r = mt.Tm_NN(reverse)

        if (gc_min <= gc_f <= gc_max and
                gc_min <= gc_r <= gc_max and
                tm_min <= tm_f <= tm_max and
                tm_min <= tm_r <= tm_max):
            return {
                "forward_primer": str(forward),
                "reverse_primer": str(reverse),
                "forward_gc": gc_f,
                "reverse_gc": gc_r,
                "forward_tm": tm_f,
                "reverse_tm": tm_r
            }
    return None


def main():
    dna_sequence = input("DNA dizisini girin: ")
    results = design_primers(dna_sequence)

    if results:
        print("İleri Primer:", results["forward_primer"])
        print("GC İçeriği:", results["forward_gc"])
        print("Erime Sıcaklığı:", results["forward_tm"])
        print("Ters Primer:", results["reverse_primer"])
        print("GC İçeriği:", results["reverse_gc"])
        print("Erime Sıcaklığı:", results["reverse_tm"])
    else:
        print("Uygun primer bulunamadı.")


if __name__ == "__main__":
    main()



