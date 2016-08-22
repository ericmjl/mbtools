from mbtools.pcr import compute_amplicon
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from hypothesis.strategies import text, integers
from hypothesis import given, assume
from entropy import shannon_entropy
import numpy as np


def test_amplify_without_overhang():
    template = seq_generator(500)
    fw_primer = Seq(template[15:30])
    re_primer = Seq(template[-30:-15]).reverse_complement()
    template = Seq(template)
    amplicon = compute_amplicon(fw_primer, re_primer, template)

    assert fw_primer in amplicon
    assert re_primer.reverse_complement() in amplicon
    assert len(amplicon) <= len(template)


def test_amplify_with_overhang():
    template = seq_generator(200)
    primer_length = 60
    ohl = 15
    fw_oh = seq_generator(ohl)
    re_oh = seq_generator(ohl)

    fw_primer = Seq(fw_oh + template[0:primer_length-ohl])
    re_primer = Seq(template[-primer_length+ohl:] + re_oh).reverse_complement()

    amplicon = compute_amplicon(fw_primer, re_primer, template)

    # This could be made better
    assert len(amplicon) == len(template) + len(fw_oh) + len(re_oh)
