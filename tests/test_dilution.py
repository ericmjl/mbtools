from mbtools.dilution import target_vol
from hypothesis import given, assume
from hypothesis.strategies import floats
from numpy import isfinite

@given(floats(), floats(), floats())
def test_target_vol(start_conc, start_vol, target_conc):
    assume(start_conc > 0)
    assume(start_vol > 0)
    assume(target_conc > 0)
    assume(isfinite(start_conc))
    assume(isfinite(start_vol))
    assume(isfinite(target_conc))

    v2 = target_vol(start_conc, start_vol, target_conc)

    assert v2 == start_conc / target_conc * start_vol
