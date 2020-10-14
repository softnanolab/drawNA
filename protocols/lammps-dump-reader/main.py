#!/usr/bin/env python
"""Manual input of LAMMPS dump and data files"""
import numpy as np
import pandas as pd
from drawNA.oxdna import System, Strand, Nucleotide
from drawNA.oxdna.system import (
    CONFIGURATION_COLUMNS, 
    TOPOLOGY_COLUMNS,
    LMP_COL_ATOMS, 
    LMP_COL_ELLIPSOIDS, 
    LMP_COL_VELOCITIES
)
from drawNA.oxdna.nucleotide import LMP_BASE, LMP_INERTIA, LMP_MASS, LMP_SHAPE

class LMPNucleotide(Nucleotide):
    def __init__(
        self, 
        _id,
        _type,
        _pos,
        _mol,
        _flag,
        _density,
        _v,
        _L,
        _shape,
        _quaternion
    ):
        _lmp_orientation = self._lmp_orientation(_quaternion)
        super().__init__(
            LMP_BASE(_type),
            _pos,
            _lmp_orientation[0],
            _lmp_orientation[2],
            _v * LMP_MASS,
            _L / LMP_INERTIA,
        )
        self._strand_index = _mol - 1
        self.index = _id - 1

    @staticmethod
    def _lmp_orientation(q) -> np.ndarray:
        """
        From pyquaternion - returns the orthonormalised
        vectors of the euclidean basis from a quaternion
        """
        # normalise first
        q = q/np.linalg.norm(q)

        quaternion_matrix = np.array([
            [q[0], -q[1], -q[2], -q[3]],
            [q[1],  q[0], -q[3],  q[2]],
            [q[2],  q[3],  q[0], -q[1]],
            [q[3], -q[2],  q[1],  q[0]]
        ])

        quaternion_bar_matrix = np.array([
            [q[0], -q[1], -q[2], -q[3]],
            [q[1],  q[0],  q[3], -q[2]],
            [q[2], -q[3],  q[0],  q[1]],
            [q[3],  q[2], -q[1],  q[0]]
        ])

        product_matrix = np.dot(
            quaternion_matrix,
            quaternion_bar_matrix.conj().transpose()
        )

        return product_matrix([1:][:, 1:])

def import_LAMMPS_data(fname: str) -> dict:
    n_atoms = int()
    n_bonds = int()
    n_strands = int()

    box = np.zeros((3, 2))
    time = int()

    lines = {}

    with open(fname, 'r') as f:
        for i, line in enumerate(f.readlines()):
            if line == 'Velocities':

    # import atom table
    atoms = pd.read_csv(fname, nrows=n_atoms, delim_whitespace=True)

    # import velocities
    velocities = pd.read_csv(fname, nrows=n_atoms, delim_whitespace=True)

    # import ellipsoid table
    ellipsoids = pd.read_csv(fname, nrows=n_atoms, delim_whitespace=True)

    # import bonds
    bonds = pd.read_csv(fname, nrows=n_bonds, delim_whitespace=True)

    dataframe = pd.DataFrame()
    assert set(dataframe.columns) == set(LMP_COL_ATOMS + LMP_COL_ELLIPSOIDS + LMP_COL_VELOCITIES)
    return {
        'n_atoms' : n_atoms,
        'n_bonds' : n_bonds,
        'n_strands' n_strands,
        'box' : box,
        'time' : time,
        'bonds' : bonds,
    }, dataframe

def detect_filetype(fname: str) -> str:
    filetype = ""
    return filetype

class LMPSystem(System):
    def __init__(self, dataframe, **kwargs):
        super().__init__(kwargs['box'], **kwargs)
        self._dataframe = dataframe
        self._meta = kwargs
        self.create_strands()
        self._initial_state = self.copy()

    def create_strands(self):
        # find n_strands
        # for each strand
        # create a new dataframe with all bonds from that strand
        # see if it is circular or not:
        ## if atom_1 has a maximum of one lonely atom
        ## and if atom_2 has a maximum of one lonely atom
        return

    def update_from_dump(self, dump : str) -> "LMPSystem":
        """Adjusts the vectors of each atom from a dump file
        and returns a new system
        """
        return
    
    def dump_to_oxDNA(self, dump_fnames : List[str]):
        """Takes a list of LAMMPS dump files and combines them
        into a single oxDNA trajectory and topology
        """
        return

    def reset(self):
        """Returns to initial configuration from LAMMPS datafile
        """
        return

    def copy(self) -> "LMPSystem":
        return

def main(fname: str, kind: str = None):
    if kind == None:
        kind = detect_filetype(fname)
    if kind == 'data':
        meta, dataframe = import_LAMMPS_data(fname)
    elif kind == 'dump':
        meta, dataframe = import_LAMMPS_dump(fname)
    system = LMPSystemFromDataFrame(dataframe, **meta)
    return

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("fname", type=str)
    parser.add_argument("")
    main(**vars(parser.parse_args()))