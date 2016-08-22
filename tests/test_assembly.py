from hypothesis.strategies import text
from hypothesis import given, assume
from hypothesis.strategies import lists, text
from mbtools.assembly import GibsonAssembler
from entropy import shannon_entropy
from .testutils import seq_generator

# @given(lists(seq_generator(500),
#              min_size=2,
#              average_size=3))


def test_gibson_assembly_class():
    """
    Most of the tests here are "data integrity" tests. The structure of the
    GibsonAssembler class has to be done right.
    """
    parts = [seq_generator(500) for i in range(3)]

    g = GibsonAssembler(parts)
    assume(len(set(parts)) > 1)  # make sure no duplicates exist
    for part in parts:
        assume(shannon_entropy(part) > 0.24)

    primer_names = ['fw_gibson', 're_gibson', '3p_sequencing', '5p_sequencing']
    for n, d in g.nodes(data=True):
        assert len(set(d.keys()).intersection(primer_names)) == 4

    assert len(g.nodes()) == len(g.edges())
    assert len(g.nodes()) == len(g.sequences)

    p = g.primers()
    assert len(p) == len(g.nodes())
    for part, primers in p.items():
        assert len(primers) == 4
