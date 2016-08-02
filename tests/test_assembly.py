from hypothesis.strategies import text
from hypothesis import given, assume
from hypothesis.strategies import lists, text
from mbtools.assembly import GibsonAssembler
from entropy import shannon_entropy


@given(lists(text(alphabet=['A', 'T', 'G', 'C'], min_size=200, max_size=500),
             min_size=2,
             average_size=3))
def test_gibson_assembly_class(parts):
    """
    Most of the tests here are "data integrity" tests. The structure of the
    GibsonAssembler class has to be done right.
    """
    g = GibsonAssembler(parts)
    assume(len(set(parts)) > 1)  # make sure no duplicates exist
    for part in parts:
        assume(shannon_entropy(part) > 0.24)

    for n, d in g.nodes(data=True):
        assert set(d.keys()) == set(['fw_gibson', 're_gibson',
                                     '3p_sequencing', '5p_sequencing'])

    assert len(g.nodes()) == len(g.edges())
    assert len(g.nodes()) == len(g.sequences)
