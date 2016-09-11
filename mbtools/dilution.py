def target_vol(c1, v1, c2):
    return c1 / c2 * v1


def input_volume(input_mass, input_conc):
    """
    Computes the required amount of volume for a given mass and concentration.

    Silently assumes that the units are:
    - input_mass: ng (nanograms)
    - input_conc: ng/µL (nanograms per microlitre)
    """
    return input_mass / input_conc


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


import collections
Sequence = collections.namedtuple("Sequence", "name dna_conc seq_length")

def cpec_equimolarity(list_seq):
    """
    Takes in PCR product lengths and PCR product concentrations,
    and returns the relative volumes of PCR product that is added to a CPEC reaction.

    Parameters:
    - list_seq: (list) a list of PCR products and their attributes in namedtuples:
        - name: (str) name given to a particular PCR product of interest
        - dna_conc: (float) in units g/µL
        - seq_length: (int) number of base **pairs**.
    Example namedtuple:
        A1 = Sequence(name="A1", dna_conc=42.3*10**-9, seq_length=1600)

    Returns:
    - rel_vol: (dict) a dictionary keyed by part name. The numbers are relative to the most concentrated piece.

    """
    d_fragment_conc = {} #dictionary of calculated concentrations of DNA sequences per ul
    relvol = {}
    for seq in list_seq:
        #Calculate the concentration of fragments of a particular sequence in solution
        fragment_conc = (float(seq.dna_conc)/(1.62*10**-21))/float(seq.seq_length)
        d_fragment_conc[str(seq.name)] = fragment_conc

    for seq_name, fragment_conc in d_fragment_conc.items():
        volume = float(max(d_fragment_conc.values()))/float(fragment_conc)
        relvol[seq_name] = volume

    return relvol
