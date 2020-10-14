#!/usr/bin/env python
"""
Script to generate an oxDNA system containing one piece of
DNA that has single-strand and double-strand block co-
polymerisation.
"""
from argparse import ArgumentParser

import numpy as np

from drawNA.oxdna import Strand, Nucleotide
from system import SystemWithFolder
from drawNA.oxdna.strand import generate_helix, FENE_EPS

def generate_system(
        n : int,
        double_strand_length : int,
        single_strand_length : int
    ) -> SystemWithFolder:
    """
    Generates an oxDNA system containing a single piece of DNA
    which has blocks of equal size of double-stranded portions
    and single-stranded portions.
    

    Parameters:
        n - number of nucleotides e.g. 100
        fraction_ds - fraction of double-strand DNA e.g. 0.5
    
    Returns:
        system - oxDNA system
    """
    # initialise system to comply with minimum image convention
    box = 2. * n * FENE_EPS
    print(f'Creating simulation system with box size: {box}')
    system = SystemWithFolder(np.array([box, box, box]))

    # create a list of strands, our main strand which we will use
    # and our complementary strand
    strands = generate_helix(n, double=True)
    
    # make a copy of the complementary strand and remove it 
    # from the main list of strands
    second_strand = strands[1].copy()
    strands = strands[:1]

    # calculate how many portions there should be
    portion_size = double_strand_length + single_strand_length
    portions = n // portion_size
    print(f'Creating {portions} portions of total size: {portion_size}')

    # iterate over all portions, adding a new complementary
    # strand to the list of strands we will use
    for i in range(portions):
        # take the nucleotides from each portion and
        # create a new strand for each double-stranded part
        # of the portion
        start = portion_size * i
        end = start + double_strand_length
        nucleotides = second_strand.nucleotides[start:end]
        new_strand = Strand(nucleotides=nucleotides)
        strands.append(new_strand)

    # finally add the strands to the system and return
    # the system
    system.add_strands(strands)
    return system

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('-n', '--number', type=int, required=True)
    parser.add_argument('-ds', '--double-stranded', type=int, required=True)
    parser.add_argument('-ss', '--single-stranded', type=int, required=True)
    parser.add_argument('-f', '--output-prefix')
    args = parser.parse_args()
    system = generate_system(
        args.number,
        args.double_stranded,
        args.single_stranded
    )
    if args.output_prefix:
        system.write_oxDNA(args.output_prefix)
