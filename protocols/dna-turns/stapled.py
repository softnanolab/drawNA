import numpy as np

from drawNA.oxdna.strand import generate_helix, POS_BACK, Strand
from drawNA.oxdna import System

FENE_LENGTH = 0.76

def main(length=17, n_strands=10, output_format='oxdna'):
    # generate a strand
    strands = []
    staples = []
    strand, staple = generate_helix(n=length, double=True)
    staples.append(staple.copy())
    strands.append(strand.copy())

    for i in range(n_strands-1):

        last_nuc = strands[-1].nucleotides[-1]
        direction = -last_nuc._a3
        a1 = -last_nuc._a1

        # ensure the backbone position is FENE_LENGTH away from
        # the backbone position of the previous nucleotide
        start = last_nuc.pos_back + (FENE_LENGTH - POS_BACK) * a1

        # generate strand above that's going in opposite direction
        strand, staple = generate_helix(
            n=length,
            start_position=start,
            direction=direction,
            a1=a1,
            double=True
        )
        strands.append(strand.copy())
        staples.append(staple.copy())

    print(staples)

    # using the two previously created strands create a new strand that
    # we will add to the system
    nucleotides = []
    for strand in strands:
        nucleotides += strand.nucleotides
    strand = Strand(nucleotides=nucleotides)

    # create system and add the final completed strand
    system = System(np.array([50., 50., 50.]))
    system.add_strand(strand)

    # create staples from staple list
    # iterate over every 2 staples in reverse
    # and append the last onto the first, and add to the system
    staples = staples[::-1]
    completed_staples = []
    for i, staple in enumerate(staples):
        if i % 2 != 0:
            continue
        try:
            new_nucleotides = []
            new_nucleotides += staples[i+1].nucleotides
            new_nucleotides += staple.nucleotides
            new_staple = Strand(nucleotides=new_nucleotides)
            system.add_strand(new_staple.copy())
        except IndexError:
            pass
    if output_format.lower() == 'oxdna':
        system.write_oxDNA('stapled_turns')
    elif output_format.lower() == 'lammps':
        system.write_lammps_data('stapled_turns')

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-l', '--length', default=16, required=False, type=int)
    parser.add_argument('-n', '--n-strands', default=10, required=False, type=int)
    parser.add_argument('-f', '--output-format', default='oxdna', required=False)
    main(**vars(parser.parse_args()))
