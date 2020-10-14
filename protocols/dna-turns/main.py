import numpy as np

from drawNA.oxdna.strand import generate_helix, POS_BACK, Strand
from drawNA.oxdna import System

FENE_LENGTH = 0.76

def main(length=16, n_strands=10):
    # generate a strand
    strands = []
    doubles = []
    strand, double = generate_helix(n=length, double=True)
    strands.append(strand.copy())
    doubles.append(double.copy())

    for i in range(n_strands-1):

        last_nuc = strands[-1].nucleotides[-1]
        direction = -last_nuc._a3
        a1 = -last_nuc._a1

        # ensure the backbone position is FENE_LENGTH away from
        # the backbone position of the previous nucleotide
        start = last_nuc.pos_back + (FENE_LENGTH - POS_BACK) * a1

        # generate strand above that's going in opposite direction
        strand, double = generate_helix(
            n=length,
            start_position=start,
            direction=direction,
            a1=a1,
            double=True,
        )
        strands.append(strand)
        doubles.append(double)

    # using the two previously created strands create a new strand that
    # we will add to the system
    nucleotides = []
    for strand in strands:
        nucleotides += strand.nucleotides
    strand = Strand(nucleotides=nucleotides)

    # create system and add the final completed strand
    system = System(np.array([20., 20., 20.]))
    system.add_strand(strand)

    actual_doubles = []
    for strand in doubles:
        nucleotides = strand.nucleotides[:5]
        actual_doubles.append(Strand(nucleotides=nucleotides))

    system.add_strands(actual_doubles)
    system.write_oxDNA('turns')

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-l', '--length', default=16, required=False, type=int)
    parser.add_argument('-n', '--n-strands', default=10, required=False, type=int)
    main(**vars(parser.parse_args()))