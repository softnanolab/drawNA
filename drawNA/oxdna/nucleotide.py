"""
Management and Storage of oxDNA nucleotides
"""
import numpy as np
import pandas as pd

from typing import List
import warnings

from .utils import get_rotation_matrix

# Emperical oxDNA constants
POS_BACK = -0.4
POS_STACK = 0.34
POS_BASE = 0.4

# Debesh's oxDNA constants
SHIFT_BASE = 0.5
SHIFT_ACROSS = 0.56 - (SHIFT_BASE * 0.9)
SHIFT_ROUND = -0.105

LMP_BASE = {
    'A' : 1,
    'C' : 2,
    'G' : 3,
    'T' : 4
}
LMP_MASS = 3.1575
LMP_INERTIA = 0.435179
LMP_SHAPE = 1.1739845031423408

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

BASE_BASE = 0.3897628551303122

class Nucleotide:
    """
    A Nucleotide is a single oxDNA particle that forms a DNA strand
    Can be added to a strand using StrandObject.add_nucleotide

    Parameters:
        base - 'A', 'T', 'C' or 'G'
        pos_com - Center of mass position vector
        a1 - Unit vector indicating orientation of backbone with respect to base
        a3 - Unit vector indicating orientation (tilting) of base with respect to backbone
        v - Linear velocity vector
        L - Angular velocity vector

    Attributes/Properties:
        _base - 'A', 'T', 'C' or 'G'
        _a1 - Unit vector indicating orientation of backbone with respect to base
        _a3 - Unit vector indicating orientation (tilting) of base with respect to backbone
        _v - Linear velocity vector
        _L - Angular velocity vector
        pos_com - centre of mass
        pos_back - backbone position
        pos_base - position of base
        pos_stack - stacking direction
        series - table entry for writing conf and top files

    """

    def __init__(
        self,
        base: str,
        pos_com: np.ndarray,
        a1: np.ndarray,
        a3: np.ndarray,
        v: np.ndarray = np.array([0.0, 0.0, 0.0]),
        L: np.ndarray = np.array([0.0, 0.0, 0.0]),
    ):
        self._base = base
        self.pos_com = pos_com
        self._a1 = a1
        self._a3 = a3
        self._v = v
        self._L = L
        self._across = None

        # make sure that the a1 and a3 vectors are orthogonal
        assert abs(round(np.dot(self._a1, self._a3))) == 0.0

        # these are accessed when the nucleotide is added
        # to an oxdna.Strand._nucleotides object
        self._strand_index = 1
        self._before = -1
        self._after = -1
        self.index = -1

    def __repr__(self) -> str:
        return f"Nucleotide[{self._base}]"

    @property
    def pos_base(self):
        """
        Returns the position of the base site
        """
        return self.pos_com + self._a1 * POS_BASE

    @property
    def pos_stack(self):
        """
        Returns the position the stacking site
        """
        return self.pos_com + self._a1 * POS_STACK

    @property
    def pos_back(self):
        """
        Returns the position of the backbone site
        """
        return self.pos_com + self._a1 * POS_BACK

    @property  # although this wasn't a property before, not sure why?
    def pos_back_rel(self):
        """
        Returns the position of the backbone centroid relative to the centre of mass
        i.e. it will be a vector pointing from the c.o.m. to the backbone
        """
        return self.pos_back - self.pos_com

    @property
    def _a2(self):
        return np.cross(self._a3, self._a1)

    @property
    def quaternion(self) -> np.ndarray:
        """
        Returns the quaternion used to create a LAMMPS configuration

        Taken from https://github.com/lorenzo-rovigatti/tacoxDNA/src/libs/oxDNA_LAMMPS.py
        """
        result = [0., 0., 0., 0.]

        a1, a2, a3 = (self._a1, self._a2, self._a3)

        q0sq = 0.25 * (a1[0] + a2[1] + a3[2] + 1.0)
        q1sq = q0sq - 0.5 * (a2[1] + a3[2])
        q2sq = q0sq - 0.5 * (a1[0] + a3[2])
        q3sq = q0sq - 0.5 * (a1[0] + a2[1])

        if q0sq >= 0.25:
            result[0] = np.sqrt(q0sq)
            result[1] = (a2[2] - a3[1]) / (4.0 * result[0])
            result[2] = (a3[0] - a1[2]) / (4.0 * result[0])
            result[3] = (a1[1] - a2[0]) / (4.0 * result[0])
        elif q1sq >= 0.25:
            result[1] = np.sqrt(q1sq)
            result[0] = (a2[2] - a3[1]) / (4.0 * result[1])
            result[2] = (a2[0] + a1[1]) / (4.0 * result[1])
            result[3] = (a1[2] + a3[0]) / (4.0 * result[1])
        elif q2sq >= 0.25:
            result[2] = np.sqrt(q2sq)
            result[0] = (a3[0] - a1[2]) / (4.0 * result[2])
            result[1] = (a2[0] + a1[1]) / (4.0 * result[2])
            result[3] = (a3[1] + a2[2]) / (4.0 * result[2])
        elif q3sq >= 0.25:
            result[3] = np.sqrt(q3sq)
            result[0] = (a1[1] - a2[0]) / (4.0 * result[3])
            result[1] = (a3[0] + a1[2]) / (4.0 * result[3])
            result[2] = (a3[1] + a2[2]) / (4.0 * result[3])

        norm = 1.0 / np.sqrt(result[0] * result[0] + result[1] * result[1] + \
                result[2] * result[2] + result[3] * result[3])

        return np.array(result) * norm

    @property
    def lammps(self) -> pd.Series:
        """
        Returns a list of Series containing the following
        information needed to write a LAMMPS configuration
        data file
        """
        result = pd.Series({
            'id' : self.index + 1,
            'type' : LMP_BASE[self._base.upper()],
            'position' : self.pos_com,
            "molecule" : self._strand_index,
            'flag' : 1,
            'density' : 1.,
            'v' : self._v / LMP_MASS,
            'L' : self._L * LMP_INERTIA,
            'shape' : [LMP_SHAPE] * 3,
            'quaternion' : self.quaternion
        })

        return result

    @property
    def series(self) -> pd.Series:
        """
        Writes a pd.Series object containing the information
        needed for writing a row in pd.DataFrame that will
        be used for writing to file.
        """
        return pd.Series(
            {
                "base": self._base,
                "position": self.pos_com,
                "a1": self._a1,
                "a3": self._a3,
                "v": self._v,
                "L": self._L,
                "strand": self._strand_index,
                "3p": self._before,
                "5p": self._after,
                "across": self._across.index,
            }
        )

    @property
    def across(self) -> int:
        if self._across:
            return self._across.index
        else:
            return -1

    @across.setter
    def across(self, nucleotide : "Nucleotide"):
        if nucleotide == None:
            self._across = None
            return

        if not isinstance(nucleotide, type(self)):
            raise TypeError(f"Set Nucleotide.across with a Nucleotide instance, you're "
                  f"using a {type(nucleotide)} instance")
        if self._across:
            warnings.warn(
                f"Setting across {nucleotide} with index {nucleotide.index}"
                f" when this nucleotide (index={self.index}) already has an"
                f" across nucleotide: {self._across} (index={self._across.index})!"
            )
            self._across.across = None
        self._across = nucleotide

    def make_3p(base: str) -> "Nucleotide":
        raise NotImplementedError

    def make_5p(
        self, 
        base: str, 
        angle: float = 0.626, 
        rise: float = 0.390,
    ) -> "Nucleotide":
        """
        Returns a new Nucleotide in the 5' direction. The angle used
        for rotating the a1 vector can be specified in radians.

        Args:
        
        base
            standard string base
        angle (0.626) 
            angle of a1 rotation around a3 axis in radians
        rise (0.390)
            difference in a3 direction between nucleotides

        Returns:
        
        A new Nucleotide generated to be in the preferred orientation
        in the 5' direction
        """
        # shift round in a1 direction
        rotation_matrix = get_rotation_matrix(self._a3, angle)
        a1 = np.dot(
                rotation_matrix,
                self._a1,
            )
        # shift up in a3 direction
        new_base = self.pos_base + rise * self._a3
        new_pos = new_base - a1 * SHIFT_BASE
        new_pos += SHIFT_ROUND * np.cross(self._a3, a1)
        return Nucleotide(base, new_pos, a1, self._a3.copy())


    def make_across(self) -> "Nucleotide":
        """
        Returns:
        
        A Nucleotide across (aka complementary)
        """
        if self._across:
            warnings.warn(
                f"Using Nucleotide.make_across when this nucleotide "
                f"(index={self.index}) already has an across nucleotide: "
                f"{self._across} (index={self._across.index})!"
            )
            self._across.across = None
        a1 = -self._a1
        a3 = -self._a3
        pos_com = self.pos_com - a1 * 2 * (SHIFT_BASE + SHIFT_ACROSS)
        nucleotide = Nucleotide(ACROSS[self._base], pos_com, a1, a3)
        nucleotide.across = self
        self.across = nucleotide
        return nucleotide


