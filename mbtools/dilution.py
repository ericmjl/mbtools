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
