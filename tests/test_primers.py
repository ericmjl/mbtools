from mbtools.primers import (anneal_temp,
                             gibson_assembly_primers,
                             gibson_sequencing_primers)
from Bio.SeqUtils.MeltingTemp import Tm_GC
from Bio.Seq import Seq
from hypothesis import given, assume
from hypothesis.strategies import floats, text
from entropy import shannon_entropy
from random import choice
from .testutils import seq_generator
import numpy as np


def test_anneal_temp():
    sequence = seq_generator(26)
    tm = Tm_GC(sequence)
    ta = anneal_temp(sequence)
    assert ta == tm - 3


@given(floats(min_value=0.01, max_value=0.1))
def test_anneal_temp_salt_corr(Na):
    sequence = seq_generator(26)
    tm = Tm_GC(sequence)
    ta = anneal_temp(sequence, salt_corr={'Na': Na})

    assert ta < tm


@given(floats(min_value=1, max_value=20))
def test_anneal_temp_chem_corr(DMSO):
    sequence = seq_generator(26)
    tm = Tm_GC(sequence)
    ta = anneal_temp(sequence, chem_corr={'DMSO': DMSO})

    assert ta < tm


def test_gibson_assembly_primers():
    part1 = seq_generator(500)
    part2 = seq_generator(300)

    primer_len = 40
    part2_fw, part1_re = gibson_assembly_primers(part1, part2,
                                                 primer_length=primer_len)

    ohl = int(round(primer_len * 1 / 3))
    ahl = primer_len - ohl

    assert part2_fw[0:ohl] == Seq(part1)[-ohl:]
    assert part2_fw[ohl:] in Seq(part2)
    assert part1_re.reverse_complement()[-ohl:] == Seq(part2)[0:ohl]
    assert part1_re.reverse_complement()[0:ahl] == Seq(part1)[-ahl:]


def test_gibson_sequencing_primers():
    part1 = seq_generator(300)
    part2 = seq_generator(900)

    up_3p, dn_5p = gibson_sequencing_primers(part1, part2)
    assert str(up_3p) in part1
    assert str(dn_5p) in part2
