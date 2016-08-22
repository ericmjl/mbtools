from Bio.SeqUtils.MeltingTemp import Tm_GC, salt_correction, chem_correction
from Bio.Seq import Seq
from Bio.Alphabet import DNAAlphabet
import numpy as np


def anneal_temp(sequence, salt_corr=None, chem_corr=None):
    """
    Annealing temperature is defined as melting temperature - 3

    salt_corr: None or dictionary of kwargs
    - Na=0, K=0, Tris=0, Mg=0, dNTPs=0, method=1, seq=None
    chem_corr: None or dictionary of kwargs
    - melting_temp, DMSO=0, fmd=0, DMSOfactor=0.75, fmdfactor=0.65,
      fmdmethod=1, GC=None
    (note: for kwargs, see Bio.SeqUtils.MeltingTemp module)
    """

    tm = Tm_GC(sequence)
    correction = 0
    if salt_corr:
        assert isinstance(salt_corr, dict),\
            "salt_corr must be a dictionary."
        correction += salt_correction(**salt_corr)

    if chem_corr:
        assert isinstance(chem_corr, dict), \
            "chem_corr must be a dictionary."
        tm = chem_correction(tm, **chem_corr)

    return tm + correction - 3


def gibson_assembly_primers(up, dn, primer_length=60):
    """
    Computes the necessary Gibson assembly primers for the junction between
    upstream (up) and downstream (dn) parts.

    `up` and `dn` should be BioPython Seq objects
    """

    assert isinstance(up, Seq)
    assert isinstance(dn, Seq)

    anl = int(round(primer_length * 2/3))  # anneal length
    ohl = int(round(primer_length * 1/3))  # overhang length

    # Downstream part's fw-assembly primer
    dn_fw = up[-ohl:] + dn[0:anl]
    # print(up[-ohl:])
    assert dn_fw[-anl:] == dn[0:anl]
    assert dn_fw[0:ohl] == up[-ohl:]
    # Upstream part's re-assembly primer
    up_re = (up[-anl:] + dn[0:ohl]).reverse_complement()
    assert up_re[0:ohl].reverse_complement() == dn[0:ohl]
    assert up_re[-anl:].reverse_complement() == up[-anl:]

    return dn_fw, up_re


def gibson_sequencing_primers(up, dn, primer_length=20):
    """
    Computes the sequencing primers necessary for verifying that junctions are
    assembled correctly.
    """
    # Upstream part's 3' junction sequencing primer
    up_3p = up[-100 - primer_length: - 100]
    assert up_3p in up
    # Downstream part's 5' junction sequencing primer.
    dn_5p = dn[100: 100 + primer_length].reverse_complement()
    assert dn_5p.reverse_complement() in dn

    return up_3p, dn_5p
