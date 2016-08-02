from mbtools.pcr import compute_amplicon
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from hypothesis.strategies import text, integers
from hypothesis import given, assume
from entropy import shannon_entropy
import numpy as np


@given(text(alphabet=['A', 'T', 'G', 'C'], min_size=200, max_size=500))
def test_amplify_without_overhang(template):
    assume(shannon_entropy(template) > 0.24)
    fw_primer = Seq(template[15:30])
    re_primer = Seq(template[-30:-15]).reverse_complement()
    template = Seq(template)
    amplicon = compute_amplicon(fw_primer, re_primer, template)

    assert fw_primer in amplicon
    assert re_primer.reverse_complement() in amplicon
    assert len(amplicon) <= len(template)


def test_amplify_with_overhang():
    """
    This test is designed as an oracle test rather than a Hypothesis test as
    above. This is because I have yet to be able to write a Hypothesis-version
    that works.
    """
    template = 'ACTACTCGCGGCATCTAGCTCACGTGACTACGTCGCGCATGGTCGTACTGTACGTCTGA' +\
               'CTATGCATCGTGTGTATTTGCGTTGGTCGTCGGCGTCATCGTAGTATTCTCTTTAGC'
    fw_oh = 'CGATTCGG'
    re_oh = 'GCGTCTTA'

    fw_primer = Seq(fw_oh + template[0:60 - len(fw_oh)])
    re_primer = Seq(template[-60 + len(re_oh):] + re_oh).reverse_complement()
    template = Seq(template)

    amplicon = compute_amplicon(fw_primer, re_primer, template)

    # This could be made better
    assert len(amplicon) == len(template) + len(fw_oh) + len(re_oh)
