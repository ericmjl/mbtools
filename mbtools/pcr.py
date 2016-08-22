"""
Author: Eric J. Ma

Purpose: A module that defines functions for doing PCR.
"""

from Levenshtein import distance
from difflib import SequenceMatcher
from Bio.Seq import Seq

import math
import numpy as np

def alignment_indices(template, primer):
    """
    Finds the optimal alignment between template and primer.

    Inputs:
    =======
    - str1, str2:  (str)

    Returns (int1, int2), (int3, int4), where:
    - int1, int2 = start, stop on str1
    - int3, int4 = start, stop on str2

    For the DNA case, we are assuming that the 5'->3' directionality of the
    two strings are identical.
    """
    s = SequenceMatcher(None, template, primer)
    m = s.find_longest_match(0, len(template), 0, len(primer))
    return (m.a, m.a + m.size), (m.b, m.b + m.size)


def compute_amplicon(fw_primer, re_primer, template):
    """
    Given a pair of fw and re primers and a template, return the amplicon.

    Inputs:
    =======
    - fw_primer, re_primer: BioPython Seq objects
    - template:             BioPython Seq objects

    It is assumed that the reverse primer is not in the same orientation as the
    template sequence.

    Returns:
    ========
    - amplicon: BioPython SeqRecord object

    Text-based illustration:

        A   B                G  H
    5-------3                5------------3
        C   D                E  F
        5-----------------------3

    Or:
        A       B                G            H
        5-------3                5------------3
        C       D                E            F
    5---------------------------------------------------3


    The basic idea here is that we need to compute the locations of D and E on
    the template.
    """
    # Identify where the fw primer anneals to template region.
    # fp   : forward primer
    # rp   : reverse primer
    # tpl  : template
    # idxs : indices
    fp_tpl_idxs, fp_idxs = alignment_indices(template,
                                             fw_primer)
    rp_tpl_idxs, rp_idxs = alignment_indices(template,
                                             re_primer.reverse_complement())

    final_sequence = fw_primer + template[fp_tpl_idxs[1]:rp_tpl_idxs[0]] +\
        re_primer.reverse_complement()

    return final_sequence


def num_cycles(fold_amp):
    """
    Computes the necessary number of cycles. Rounds up to the next integer.
    """
    return math.ceil(np.log10(fold_amp) / np.log10(2))


def num_cycles(fold_amp):
    """
    Computes the necessary number of cycles. Rounds up to the next integer.
    """
    return math.ceil(np.log10(fold_amp) / np.log10(2))


def input_plasmid_mass(target_len, plasmid_len, target_mass):
    """
    Computes the mass of total input plasmid mass required to get a given mass
    of target DNA.

    Silently assumes that the units are:
    - target_len: bp (base pairs)
    - plasmid_len: bp (base pairs)
    - target_mass: ng (nanograms)
    """
    return target_mass * plasmid_len / target_len


def input_volume(input_mass, input_conc):
    """
    Computes the required amount of volume for a given mass and concentration.

    Silently assumes that the units are:
    - input_mass: ng (nanograms)
    - input_conc: ng/µL (nanograms per microlitre)
    """
    return input_mass / input_conc


def num_cycles(fold_amp):
    """
    Computes the necessary number of cycles. Rounds up to the next integer.
    """
    return math.ceil(np.log10(fold_amp) / np.log10(2))
