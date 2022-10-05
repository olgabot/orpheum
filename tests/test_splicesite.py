import pytest

from orpheum.splicesite import SpliceSite


@pytest.fixture
def read_with_3prime_splice_site():
    return 'AGCATATCTGTAAGCTGCATAGCCTCCCTGCAGGGGAGGTGGGAAGGAGAGATGGCAGATGGATTTAGTGAGGCTCGCTGTTGCTGAAGC'


@pytest.fixture
def read_with_3prime_splice_site_reverse_complement():
    return 'GCTTCAGCAACAGCGAGCCTCACTAAATCCATCTGCCATCTCTCCTTCCCACCTCCCCTGCAGGGAGGCTATGCAGCTTACAGATATGCT'


@pytest.fixture
def read_with_5prime_splice_site():
    return 'TCCTGTCTGCCTGGACTCACCCTCCAAAGCCCACCACAACCCCCGTGTAGCCTGCTGCTCTTAGCCTTACCTGTCACTGTAGGCTGCTGC'


@pytest.fixture
def read_with_short_exon():
    """Contains both splice sites"""
    return 'TTTTTCCCCTCTCCTACTGTCTGGATTTCTTGTTCAACAGAAATGAAATACCTAGCAGGTACTGCCATTGGGCACGGTGTGGGGTGTGCC'


@pytest.fixture
def splice_site():
    return SpliceSite()


@pytest.mark.usefixtures("seq")
def test_score3_full_seq(seq, splice_site):
    score3 = splice_site.score3_full_seq(seq)
    assert score3 == {'AAATCGCAATATAACTGTAAATC': -7.102585997971779,
                      'TGCTTAATACTGACATCAATAAT': -9.817013352622899,
                      'CATCAATAATATTAGGAAAATCG': -17.414368922924915,
                      'AATAATATTAGGAAAATCGCAAT': -24.799525308177547,
                      'GCAATATAACTGTAAATCCTGTT': -14.20507047687658,
                      'ATCGCAATATAACTGTAAATCCT': -9.610111135110257,
                      'CTGACATCAATAATATTAGGAAA': -8.128742704102592,
                      'TAATACTGACATCAATAATATTA': -18.76860799830603,
                      'AATATAACTGTAAATCCTGTTCT': -14.259589459115302,
                      'CAATAATATTAGGAAAATCGCAA': -17.195161307409496,
                      'AATATTAGGAAAATCGCAATATA': -17.93568502188441,
                      'ATACTGACATCAATAATATTAGG': -19.38583094680745,
                      'GACATCAATAATATTAGGAAAAT': -21.703495638365663,
                      'GCTTAATACTGACATCAATAATA': -18.98106450410893,
                      'TAATATTAGGAAAATCGCAATAT': -14.63751869602015,
                      'CGCTTGCTTAATACTGACATCAA': -6.7113228113808105,
                      'AAAATCGCAATATAACTGTAAAT': -24.33886291608301,
                      'GGAAAATCGCAATATAACTGTAA': -11.069537633194146,
                      'TTAGGAAAATCGCAATATAACTG': -9.76864831635461,
                      'ATAATATTAGGAAAATCGCAATA': -32.23638784350524,
                      'TCAATAATATTAGGAAAATCGCA': -23.80532020973892,
                      'GAAAATCGCAATATAACTGTAAA': -16.79865682393187,
                      'TATTAGGAAAATCGCAATATAAC': -11.596392203495267,
                      'AATCGCAATATAACTGTAAATCC': -12.310014181581165,
                      'CTTGCTTAATACTGACATCAATA': -13.453928592776515,
                      'TATAACTGTAAATCCTGTTCTGT': -11.78807238392215,
                      'ACTGACATCAATAATATTAGGAA': 2.828745710625767,
                      'CAATATAACTGTAAATCCTGTTC': -8.593974827356204,
                      'TCGCAATATAACTGTAAATCCTG': -17.93209369313127,
                      'AGGAAAATCGCAATATAACTGTA': -18.357938319845392,
                      'TTAATACTGACATCAATAATATT': -9.066921638439034,
                      'ACATCAATAATATTAGGAAAATC': -14.950883399873405,
                      'CGCAATATAACTGTAAATCCTGT': -15.50730119239035,
                      'AATACTGACATCAATAATATTAG': -8.348847125772568,
                      'TTGCTTAATACTGACATCAATAA': -10.269859396055795,
                      'TGACATCAATAATATTAGGAAAA': -20.710077578368605,
                      'CTTAATACTGACATCAATAATAT': -7.322704706892016,
                      'TAGGAAAATCGCAATATAACTGT': -12.926789916980802,
                      'ATAACTGTAAATCCTGTTCTGTC': -11.219404402763239,
                      'TACTGACATCAATAATATTAGGA': -15.398442625465771,
                      'ATTAGGAAAATCGCAATATAACT': -24.736060400137177,
                      'GCTTGCTTAATACTGACATCAAT': -18.536950834036194,
                      'ATCAATAATATTAGGAAAATCGC': -17.993479571230367,
                      'ATATAACTGTAAATCCTGTTCTG': -24.0852128280598,
                      'ATATTAGGAAAATCGCAATATAA': -23.992456230714808}


@pytest.mark.usefixtures("seq")
def test_score5_full_seq(seq, splice_site):
    score5 = splice_site.score5_full_seq(seq)
    assert score5 == {'ATAATATTA': -16.265973552227862,
                      'ATTAGGAAA': -24.820700025039677,
                      'GACATCAAT': -17.14788039505562, 'AAATCGCAA': -23.511371862754768,
                      'CTGACATCA': -14.33904322205772, 'CAATATAAC': -24.059209440299703,
                      'ATCGCAATA': -9.409699579451123, 'AACTGTAAA': -31.389703464795126,
                      'CTTAATACT': -31.89917927332821, 'AGGAAAATC': -11.743472809863412,
                      'ACTGTAAAT': 0.40684641320162673,
                      'TACTGACAT': -21.332563499848952,
                      'ATATTAGGA': -8.747683453622983, 'GAAAATCGC': -20.281754106484698,
                      'ATACTGACA': -16.43416355011715, 'TAATATTAG': -30.821585097320686,
                      'AAATCCTGT': -18.105391209707104,
                      'TGTTCTGTC': -47.489518790465674,
                      'GGAAAATCG': -22.757089933507622,
                      'CTGTTCTGT': -5.4522449557897605,
                      'CTTGCTTAA': -31.588144853194635, 'TAGGAAAAT': -1.562461279006692,
                      'TTAATACTG': -16.613607952076425, 'TATTAGGAA': -32.29059378121649,
                      'GTAAATCCT': -34.10424699253535, 'AATATAACT': -2.2684909357219722,
                      'TCAATAATA': -11.207415682661289,
                      'ATCCTGTTC': -21.461832985147243,
                      'TCGCAATAT': -12.039999479378713,
                      'TATAACTGT': -21.811591377110307,
                      'GCAATATAA': -15.03983265052683, 'TAAATCCTG': -27.037682075033217,
                      'ATAACTGTA': -36.762373450585315,
                      'TAATACTGA': -26.433483540542618,
                      'TGACATCAA': -39.06277325542777, 'TGCTTAATA': -10.190118292928808,
                      'TTAGGAAAA': -13.021019490331225,
                      'CTGTAAATC': -13.240650319057707,
                      'ACATCAATA': -18.619274083812936,
                      'ATATAACTG': -22.927341136520752,
                      'CCTGTTCTG': -23.367831801039955,
                      'AATCGCAAT': -23.040602972672993,
                      'ACTGACATC': -24.439729286412696,
                      'AAAATCGCA': -24.434788594379036,
                      'AATATTAGG': -14.459684614680382,
                      'CAATAATAT': -14.712041512276047,
                      'CGCTTGCTT': -13.136713399814877, 'AATACTGAC': -33.09686686274021,
                      'AATCCTGTT': -31.75685997477841, 'TAACTGTAA': -20.93390597921523,
                      'AATAATATT': -27.54647927552244,
                      'GCTTAATAC': -23.68828569986888, 'TGTAAATCC': -21.4958542881205,
                      'GCTTGCTTA': -38.94716476784485,
                      'TTGCTTAAT': -10.219833834412416, 'CATCAATAA': -19.86663372492116,
                      'ATCAATAAT': -30.76581982995747, 'TCCTGTTCT': -35.913430128402,
                      'CGCAATATA': -34.24815708622815}


def test_argmax_score5(seq, splice_site):
    kmer5, score5 = splice_site.argmax_score5(seq)
    assert kmer5 == 'ACTGTAAAT'
    assert score5 == 0.40684641320162673


def test_argmax_score3(seq, splice_site):
    kmer3, score3 = splice_site.argmax_score3(seq)
    assert kmer3 == 'ACTGACATCAATAATATTAGGAA'
    assert score3 == 2.828745710625767


def test_extract_exonic(read_with_3prime_splice_site, splice_site):
    exonic = splice_site.extract_exonic(read_with_3prime_splice_site)
    assert 0, exonic


def test_extract_exonic(read_with_5prime_splice_site, splice_site):
    exonic = splice_site.extract_exonic(read_with_5prime_splice_site)
    assert 0, exonic


def test_extract_exonic(read_with_short_exon, splice_site):
    exonic = splice_site.extract_exonic(read_with_short_exon)
    assert exonic == 'AAATGAAATACCTAGCAG'
