from maxentpy import maxent_fast
from orpheum.compare_kmer_content import kmerize


class SpliceSite:
    def __init__(self):
        self.maxent_matrix5 = maxent_fast.load_matrix(5)
        self.maxent_matrix3 = maxent_fast.load_matrix(3)

    def score3_single_23mer(self, three_prime_23mer):
        return maxent_fast.score3(three_prime_23mer, matrix=self.maxent_matrix3)

    def score5_single_9mer(self, five_primer_9mer):
        return maxent_fast.score5(five_primer_9mer, matrix=self.maxent_matrix5)

    def score5_full_seq(self, sequence):
        maxent_5p_scores = {kmer: self.score5_single_9mer(kmer) for kmer in kmerize(sequence, 9)}
        return maxent_5p_scores

    def score3_full_seq(self, sequence):
        maxent_3p_scores = {kmer: self.score3_single_23mer(kmer) for kmer in kmerize(sequence, 23)}
        return maxent_3p_scores
