"""
Storage and Management of oxDNA strands via the Strand
class
"""
import numpy as np
import pandas as pd
from typing import List
import re
import random
from copy import deepcopy

from .nucleotide import Nucleotide
from .utils import get_rotation_matrix, round_to_multiple

# Constants
PI = np.pi

# Emperical oxDNA constants
POS_BACK = -0.4
POS_STACK = 0.34
POS_BASE = 0.4

FENE_EPS = 2.0
FENE_LENGTH = 0.76

# Center of the double strand
CM_CENTER_DS = POS_BASE + 0.2

# Ideal distance of base sites of 2 nts (to be base paired in a duplex)
BASE_BASE = 0.3897628551303122

number_to_base = {0: "A", 1: "G", 2: "C", 3: "T"}
base_to_number = {"A": 0, "G": 1, "C": 2, "T": 3}


class Strand:
    """
    Collection of nucleotides in the 3' -> 5' direction
    Can be added to the system using SystemObject.add_strand

    Parameters:
        nucleotides ([]) - list of nucleotides to form the strand
    
    Attributes/Properties:
        index - strand_id
        sequence - string sequence of bases
        nucleotides - list of Nucleotides
        dataframe - table of values for top and conf files
        copy - get a copy of a Strand instance
    """

    def __init__(self, nucleotides: list = None):
        if nucleotides:
            self._nucleotides = nucleotides
        else:
            self._nucleotides = []

        self.index = 1
        self._nucleotide_shift = 0
        self.circular = False

    def __repr__(self):
        return f"3>[{self.index}]>{self.sequence}>>5"

    @property
    def sequence(self) -> str:
        if len(self.nucleotides) == 0:
            return ""

        return "".join([i._base for i in self.nucleotides])

    @property
    def nucleotides(self) -> list:
        """
        Returns the list of nucleotides, and allocates
        the strand index to each nucleotide for its 
        Nucleotide.series property.
        """
        if len(self._nucleotides) == 0:
            return []

        for i, nucleotide in enumerate(self._nucleotides):
            nucleotide.index = i + self._nucleotide_shift
            if i == 0:
                #nucleotide._before = -1
                if self.circular:
                    nucleotide._before = \
                        len(self._nucleotides) + self._nucleotide_shift - 1
                else:
                    nucleotide._before = -1

            else:
                nucleotide._before = i - 1 + self._nucleotide_shift
            if i == len(self._nucleotides) - 1:
                #nucleotide._after = -1
                if self.circular:
                    nucleotide._after = self._nucleotide_shift
                else:
                    nucleotide._after = -1

            else:
                nucleotide._after = i + 1 + self._nucleotide_shift
            nucleotide._strand_index = self.index

        #if self.circular:
            #self._nucleotides[-1].after = self._nucleotides[0].index
            #self._nucleotides[0].before = self._nucleotides[-1].index


        #    self._nucleotides[0].before = self._nucleotides[-1]#.index
        #    self._nucleotides[-1].after = self._nucleotides[0]#.index
        #    print('connected both ends of the strand to form a loop.')

        return self._nucleotides

    def add_nucleotide(self, nucleotide: Nucleotide):
        self._nucleotides.append(nucleotide)

    @property
    def dataframe(self) -> pd.DataFrame:
        """
        Returns a pd.Dataframe which can be used to write
        the information for both the configuration and
        topology files.
        """
        result = pd.DataFrame([i.series for i in self.nucleotides])
        return result

    @property
    def bonds(self) -> pd.DataFrame:
        """
        Returns a dataframe containing the information needed to write
        the bonds section of a LAMMPS configuration data file
        """
        result = pd.DataFrame(
            {
                "type": [1] * len(self.nucleotides),
                "atom_1": [i.index + 1 for i in self.nucleotides],
                "atom_2": [i._after + 1 for i in self.nucleotides],
            }
        )
        result = result[result.atom_1 != 0]
        result = result[result.atom_2 != 0]
        return result

    @property
    def lammps(self) -> List[pd.DataFrame]:
        """
        Returns a list of DataFrames containing the following
        information needed to write a LAMMPS configuration
        data file
        """
        result = pd.DataFrame([i.lammps for i in self.nucleotides])
        return result

    def __len__(self) -> int:
        return len(self._nucleotides)

    @sequence.setter
    def sequence(self, seq: str):
        seq = seq.upper()
        if len(seq) > len(self._nucleotides):
            raise ValueError(
                f"Sequence length must = strand length {len(self._nucleotides)}"
            )
        if re.findall("[^ACGT]", seq):
            raise ValueError("Sequence can only contain A, G, C or T")

        for i, base in enumerate(seq):
            self._nucleotides[i]._base = base

    # When a strand is added using add_strand in system.py
    # We create a new object, so the strand is not accessible
    # from anywhere other than the system
    def copy(self):
        return deepcopy(self)

    def transform(self, matrix: np.ndarray):
        for nucleotide in self._nucleotides:
            nucleotide.transform(matrix)
        return

    def translate(self, translation_vector: np.ndarray):
        for nucleotide in self._nucleotides:
            nucleotide.translate(translation_vector)
        return

    def rotate(self, rotator: np.ndarray):
        for nucleotide in self._nucleotides:
            nucleotide.rotate(rotator)
        return


def generate_helix(
    n: int = None,
    sequence: str = None,
    start_position: np.ndarray = np.array([0., 0., 0.]),
    direction: np.ndarray = np.array([0., 1., 0.]),
    a1: np.ndarray = np.array([1., 0., 0.]),
    initial_rotation: float = None,
    double: bool = False,
    enforce_180: bool = True,
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
        a1 = np.dot(get_rotation_matrix(direction, initial_rotation), a1)

    if enforce_180 and not n == 1:
        half_turns = round_to_multiple(n/bp_per_turn, 0.5, 1)
        # minus one because for x nt there will only be x-1 bonds
        angle = np.radians(360*half_turns / (n-1) )
    else:
        angle = 0.626

    # initialise strand list
    strands = []

    # create main strand
    strand = Strand()
    strand.add_nucleotide(
        Nucleotide(sequence[0], start_position, a1=a1, a3=direction))

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
