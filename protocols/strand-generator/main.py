import numpy as np
import random
import sys

from drawNA.oxdna import Strand, Nucleotide, System
from drawNA.oxdna.strand import BASE_BASE, FENE_EPS, POS_BACK, POS_STACK, generate_helix

POS_BASE = 0.5
BASE_BASE = 0.56 - (POS_BASE * 0.9)

import scipy.linalg as la

def get_rotation_matrix(axis, theta):
    """Rotation matrix around a vector"""
    return la.expm(np.cross(np.eye(3), axis/la.norm(axis)*theta))

ACROSS = {
    'A' : 'T',
    'C' : 'G',
    'G' : 'C',
    'T' : 'A',
    'a' : 't',
    't' : 'a',
    'c' : 'g',
    'g' : 'c'
}

def get_3p(nuc : Nucleotide, base : str = 'A') -> Nucleotide:
    return

def get_5p(nuc : Nucleotide, base : str = 'T') -> Nucleotide:
    # shift round in a1 direction
    angle = (np.pi / 180.0) * 35.9
    rotation_matrix = get_rotation_matrix(nuc._a3, angle)
    a1 = np.dot(
            rotation_matrix,
            nuc._a1,
        )
    # shift up in a3 direction
    new_base = nuc.pos_base + 0.39 * nuc._a3
    new_pos = new_base - a1 * (POS_BASE) - 0.105 * np.cross(nuc._a3, a1)
    return Nucleotide(base, new_pos, a1, nuc._a3.copy())

def get_across(nuc : Nucleotide) -> Nucleotide:
    """
    Returns Nucleotide o
    """
    a1 = -nuc._a1
    a3 = -nuc._a3
    pos_com = nuc.pos_com - a1 * 2 * (POS_BASE + BASE_BASE)
    return Nucleotide(ACROSS[nuc._base], pos_com, a1, a3)

def main(number : int = 10, double : bool = False):
    # Add Strand 1 - ssDNA
    n = number
    nucleotides = []
    print("Creating a nucleotide:")
    nucleotides.append(
        Nucleotide(
            random.choice(['A', 'T', 'C', 'G']),
            np.array([0., -10., 0.]),
            a1 = np.array([1., 0., 0.]),
            a3 = np.array([0., 1., 0.])
        )
    )
    print(f"Nucleotide #0: {nucleotides[0]}")
    print("Creating more nucleotides...\n")
    for i in range(n):
        nucleotides.append(get_5p(nucleotides[-1], base=random.choice(['A', 'T', 'C', 'G'])))
    strand = Strand(nucleotides=nucleotides)
    print(f"Strand: {strand}")
    system = System(np.array([20., 20., 20.]))
    system.add_strand(strand)
    # Add Strand 2 - complementary ssDNA -> dsDNA
    if double:
        nucleotides = []
        for nucleotide in strand.nucleotides[::-1]:
            nucleotides.append(get_across(nucleotide))
    #    for i in range(10):
    #        nucleotides.append(get_5p(nucleotides[-1]))
        strand = Strand(nucleotides=nucleotides)
        print(f"Adding double: {strand}")
        system.add_strand(strand)
    system.write_oxDNA('local')
    return

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-f', '--output-format', required=False, default='oxdna')
    print('Running Strand Generator...\n')
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-n', '--number', type=int, required=False, default=10)
    parser.add_argument('--double', default=False, required=False, dest='double', action='store_true')
    parser.set_defaults(double=False)
    args = vars(parser.parse_args())
    print(args)

    main(**args)
