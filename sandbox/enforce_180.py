from drawNA.oxdna import System, Strand, Nucleotide, POS_BACK, get_rotation_matrix
from drawNA.lattice.utils import round_to_multiple
import numpy as np
from typing import List
import random

FENE_LENGTH = 0.76

def generate_helix(
    n: int = None,
    sequence: str = None,
    start_position: np.ndarray = np.array([0., 0., 0.]),
    direction: np.ndarray = np.array([0., 1., 0.]),
    a1: np.ndarray = np.array([1., 0., 0.]),
    initial_rotation: float = None,
    double: bool = False,
    enforce_180: bool = False,
    bp_per_turn: float = 10.45,
) -> List[Strand]:

    # handle sequence/n arguments
    if sequence and n:
        if len(sequence) != n:
            difference = len(sequence) - n

            # sequence longer than n
            if difference > 0:
                sequence = sequence[:n]

            # n longer than sequence
            else:
                sequence += ''.join(
                    [random.choice(['A', 'T', 'C', 'G']) for i in range(-difference)]
                )

    elif sequence:
        n = len(sequence)
    elif n:
        sequence = ''.join([random.choice(['A', 'T', 'C', 'G']) for i in range(n)])
    else:
        raise TypeError(
            'Please provide either the number of base-pairs or a sequence'
        )

    # handle a1/angle arguments
    if initial_rotation:
        a1 = get_rotation_matrix(direction, initial_rotation)

    if enforce_180:
        half_turns = round_to_multiple(n/bp_per_turn, 0.5, 1)
        angle = np.radians(360*half_turns / (n-1) )
    else:
        angle = 0.626

    # initialise strand list
    strands = []

    # create main strand
    strand = Strand()
    strand.add_nucleotide(
        Nucleotide(sequence[0], start_position, a1=a1, a3=direction)
        )

    for base in sequence[1:]:
        strand.add_nucleotide(strand.nucleotides[-1].make_5p(base, angle))

    # add to strand list which will be returned
    strands.append(strand.copy())

    # create across strand 
    if double:
        strand = Strand()
        # iterate over nucleotides from original strand but in reverse
        for nucleotide in strands[0].nucleotides[::-1]:
            strand.add_nucleotide(nucleotide.make_across())
        strands.append(strand.copy())
    
    return strands

def main(length=16, n_strands=10):
    # generate a strand
    strands = []
    strand = generate_helix(n=length, enforce_180= True)[0]
    strands.append(strand.copy())

    for i in range(n_strands-1):

        last_nuc = strands[-1].nucleotides[-1]
        direction = -last_nuc._a3
        a1 = -last_nuc._a1

        # ensure the backbone position is FENE_LENGTH away from
        # the backbone position of the previous nucleotide
        start = last_nuc.pos_back + (FENE_LENGTH - POS_BACK) * a1
        # generate strand above that's going in opposite direction
        strand = generate_helix(
            n=length,
            start_position=start,
            direction=direction,
            a1=a1,
            enforce_180=True,
        )[0]
        strands.append(strand)

    # using the two previously created strands create a new strand that
    # we will add to the system
    nucleotides = []
    for strand in strands:
        nucleotides += strand.nucleotides
    strand = Strand(nucleotides=nucleotides)

    # create system and add the final completed strand
    system = System(np.array([50., 50., 50.]))
    system.add_strand(strand)
    system.write_oxDNA('turns_18nt')
    return system

if __name__ == '__main__':
    system = main(18)