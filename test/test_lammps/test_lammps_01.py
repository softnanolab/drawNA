from drawNA.oxdna import System, Strand, Nucleotide
from drawNA.oxdna.strand import generate_helix
import numpy as np
import pandas as pd


def test_nucleotide():
    print(Nucleotide(
                "A",
                np.array([1.0, 0.0, 0.0]),
                np.array([1.0, 0.0, 0.0]),
                np.array([0.0, 0.0, 1.0]),
            ).lammps)
    return

def test_strand():
    strand = Strand(
        [
            Nucleotide(
                "A",
                np.array([1.0, 0.0, 0.0]),
                np.array([1.0, 0.0, 0.0]),
                np.array([0, 0.0, 1.0]),
            ),
            Nucleotide(
                "A",
                np.array([2.0, 0.0, 0.0]),
                np.array([1.0, 0.0, 0.0]),
                np.array([0, 0.0, 1.0]),
            ),
            Nucleotide(
                "A",
                np.array([3.0, 0.0, 0.0]),
                np.array([1.0, 0.0, 0.0]),
                np.array([0, 0.0, 1.0]),
            ),
            Nucleotide(
                "A",
                np.array([4.0, 0.0, 0.0]),
                np.array([1.0, 0.0, 0.0]),
                np.array([0, 0.0, 1.0]),
            ),
            Nucleotide(
                "A",
                np.array([5.0, 0.0, 0.0]),
                np.array([1.0, 0.0, 0.0]),
                np.array([0, 0.0, 1.0]),
            ),
            Nucleotide(
                "A",
                np.array([6.0, 0.0, 0.0]),
                np.array([1.0, 0.0, 0.0]),
                np.array([0, 0.0, 1.0]),
            ),
        ]
    )
    print(strand.lammps)
    return

def test_system():
    system = System(np.array([50., 50., 50.]))
    system.add_strands(generate_helix(10, double=True))
    print(system.lammps)
    print(system.bonds)
    
    return system

if __name__ == '__main__':
    test_nucleotide()
    test_strand()
    system = test_system()
    system.write_LAMMPS()