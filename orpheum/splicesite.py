from maxentpy import maxent_fast

from orpheum.compare_kmer_content import kmerize

old_chars = "ACGT"
replace_chars = "TGCA"
REVERSE_COMPLEMENT = str.maketrans(old_chars, replace_chars)


def reverse_complement(sequence):
    return sequence.translate(REVERSE_COMPLEMENT)[::-1]


class SpliceSite:
    def __init__(self):
        self.maxent_matrix5 = maxent_fast.load_matrix(5)
        self.maxent_matrix3 = maxent_fast.load_matrix(3)

    def score3_single_23mer(self, three_prime_23mer):
        return maxent_fast.score3(three_prime_23mer, matrix=self.maxent_matrix3)

    def score5_single_9mer(self, five_primer_9mer):
        return maxent_fast.score5(five_primer_9mer, matrix=self.maxent_matrix5)

    def score5_full_seq(self, sequence):
        maxent_5p_scores = {kmer: self.score5_single_9mer(kmer) for kmer in
                            kmerize(sequence, 9)}
        return maxent_5p_scores

    def score3_full_seq(self, sequence):
        maxent_3p_scores = {kmer: self.score3_single_23mer(kmer) for kmer in
                            kmerize(sequence, 23)}
        return maxent_3p_scores

    def argmax_score5(self, sequence):
        score5 = self.score5_full_seq(sequence)
        argmax_kmer = max(score5, key=score5.get)
        return argmax_kmer, score5[argmax_kmer]

    def argmax_score3(self, sequence):
        score3 = self.score3_full_seq(sequence)
        argmax_kmer = max(score3, key=score3.get)
        return argmax_kmer, score3[argmax_kmer]

    def argmax_score5_with_rc(self, seq, seq_rc):
        kmer5, score5 = self.argmax_score5(seq)
        kmer5_rc, score5_rc = self.argmax_score5(seq_rc)
        if score5 > score5_rc:
            return kmer5, score5, 'forward'
        else:
            return kmer5_rc, score5_rc, 'reverse'

    def argmax_score3(self, seq, seq_rc):
        kmer3, score3 = self.argmax_score3(seq)
        kmer3_rc, score3_rc = self.argmax_score3(seq_rc)
        if score3 > score3_rc:
            return kmer3, score3, 'forward'
        else:
            return kmer3_rc, score3_rc, 'reverse'

    def extract_exonic(self, sequence):
        seq_rc = reverse_complement(sequence)
