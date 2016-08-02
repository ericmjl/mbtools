from mbtools.primers import anneal_temp, gibson_assembly_primers
from Bio.SeqUtils.MeltingTemp import Tm_GC
from Bio.Seq import Seq
from hypothesis import given, assume
from hypothesis.strategies import floats, text
from entropy import shannon_entropy


def test_anneal_temp():
    sequence = 'atcgctacgcggctagcgtctagcat'

    tm = Tm_GC(sequence)

    ta = anneal_temp(sequence)

    assert ta == tm - 3


@given(floats(min_value=0.01, max_value=0.1),
       text(alphabet=['A', 'T', 'G', 'C'], min_size=15, max_size=60))
def test_anneal_temp_salt_corr(Na, sequence):
    tm = Tm_GC(sequence)
    ta = anneal_temp(sequence, salt_corr={'Na': Na})

    assert ta < tm


@given(floats(min_value=1, max_value=20),
       text(alphabet=['A', 'T', 'G', 'C'], min_size=15, max_size=60))
def test_anneal_temp_chem_corr(DMSO, sequence):
    tm = Tm_GC(sequence)
    ta = anneal_temp(sequence, chem_corr={'DMSO': DMSO})

    assert ta < tm


@given(text(alphabet=['A', 'T', 'G', 'C'], min_size=200, max_size=1000),
       text(alphabet=['A', 'T', 'G', 'C'], min_size=200, max_size=1000))
def test_gibson_assembly_primers(part1, part2):
    assume(shannon_entropy(part1) > 0.24)
    assume(shannon_entropy(part2) > 0.24)

    part2_fw, part1_re = gibson_assembly_primers(part1, part2)

    assert part2_fw[15:] in Seq(part2)
    assert part1_re[15:].reverse_complement() in Seq(part1)
