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

def cpec_equimolarity(dict_conc_bp):
    """
    For every sequence, the key is the name of the sequence, and the value is a list of
    two numbers: 1) concentration of DNA in solution (g/ul) and
    2) the number of base pairs in the sequence
    """
    d_concentrations = {}
    l_concentrations = []
    relvol = {}
    for key, value in dict_conc_bp.items():
        #Calculate the concentration of fragments of a particular sequence in solution
        fragmentperul = (float(value[0])/(1.62*10**-21))/float(value[1])

        d_concentrations[key] = fragmentperul
        l_concentrations.append(fragmentperul)
    for key, value in d_concentrations.items():
        """
        By dividing the highest concentration by each of the other concentrations,
        we know how much volume of a given sequence needs to used to make it equimolar to
        one unit volume of the sequence with the highest concentration
        """

        volume = float(max(l_concentrations))/float(value)
        relvol[key] = volume

    return relvol
