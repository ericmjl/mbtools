"""
Author: Eric J. Ma

Purpose: A module that defines functions for doing PCR.
"""

from Levenshtein import distance
from difflib import SequenceMatcher
from Bio.Seq import Seq


def alignment_indices(str1, str2):
    """
    Finds the optimal alignment between two strings.

    Inputs:
    =======
    - str1, str2:  (str)

    Returns (int1, int2), (int3, int4), where:
    - int1, int2 = start, stop on str1
    - int3, int4 = start, stop on str2

    For the DNA case, we are assuming that the 5'->3' directionality of the
    two strings are identical.
    """
    s = SequenceMatcher(None, str1, str2)
    m = s.find_longest_match(0, len(str1), 0, len(str2))

    return (m.a, m.a + m.size), (m.b, m.b + m.size)


def compute_amplicon(fw_primer, re_primer, template):
    """
    Given a pair of fw and re primers and a template, return the amplicon.

    Inputs:
    =======
    - fw_primer, re_primer: BioPython SeqRecord objects
    - template: BioPython SeqRecord objects

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
    fp_idxs, fp_tpl_idxs = alignment_indices(fw_primer,
                                             template)
    rp_idxs, rp_tpl_idxs = alignment_indices(re_primer.reverse_complement(),
                                             template)

    final_sequence = \
        fw_primer +\
        template[fp_tpl_idxs[1]:rp_tpl_idxs[0]] +\
        re_primer.reverse_complement()

    return final_sequence
