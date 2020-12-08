from drawNA.oxdna import Strand, Nucleotide, System
from drawNA.oxdna.strand import generate_helix, POS_BACK

FENE_LENGTH = 0.76

#Generation of initial strand configuration
def generate_system(box, length=16, n_strands=10, stapled=5):
    strands = []
    doubles = []
    strand, double = generate_helix(n=length, double=True)
    strands.append(strand.copy())
    doubles.append(double.copy())

    for i in range(n_strands-1):

        last_nuc = strands[-1].nucleotides[-1]
        direction = -last_nuc._a3
        a1 = -last_nuc._a1

        # ensure the backbone position is FENE_LENGTH away from the backbone position of the previous nucleotide
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

        # using the two previously created strands to create a new strand that is added to the system
        nucleotides = []
        for strand in strands:
            nucleotides += strand.nucleotides

        strand = Strand(nucleotides=nucleotides)

        # create system and add the final completed strand
        main_system = System(box)
        main_system.add_strand(strand)

        actual_doubles = []
        for strand in doubles:
            nucleotides = strand.nucleotides[:stapled]
            actual_doubles.append(Strand(nucleotides=nucleotides))

        main_system.add_strands(actual_doubles)

    return main_system