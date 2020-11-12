from drawNA.oxdna import System, Strand, Nucleotide
import numpy as np
import pandas as pd


def test_System():
    system = System(np.array([50.0, 50.0, 50.0]))
    strand_1 = Strand(
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
    strand_2 = Strand(
        [
            Nucleotide(
                "T",
                np.array([1.0, 2.0, 0.0]),
                np.array([1.0, 0.0, 0.0]),
                np.array([0, 0.0, 1.0]),
            ),
            Nucleotide(
                "T",
                np.array([2.0, 2.0, 0.0]),
                np.array([1.0, 0.0, 0.0]),
                np.array([0, 0.0, 1.0]),
            ),
            Nucleotide(
                "T",
                np.array([3.0, 2.0, 0.0]),
                np.array([1.0, 0.0, 0.0]),
                np.array([0, 0.0, 1.0]),
            ),
            Nucleotide(
                "T",
                np.array([4.0, 2.0, 0.0]),
                np.array([1.0, 0.0, 0.0]),
                np.array([0, 0.0, 1.0]),
            ),
            Nucleotide(
                "T",
                np.array([5.0, 2.0, 0.0]),
                np.array([1.0, 0.0, 0.0]),
                np.array([0, 0.0, 1.0]),
            ),
            Nucleotide(
                "T",
                np.array([6.0, 2.0, 0.0]),
                np.array([1.0, 0.0, 0.0]),
                np.array([0, 0.0, 1.0]),
            ),
        ]
    )
    system.add_strands([strand_1, strand_2])

    assert isinstance(system, System)

    assert system.E_tot == 0.0
    assert len(system.dataframe.count()) == 10
    assert len(system.configuration.count()) == 5
    assert len(system.topology.count()) == 4

    assert len(system.strands) == 2
    assert len(system.nucleotides) == 12

    system.strands[0].sequence = "AGAGAG"
    system.add_strand(system.strands[0].copy())
    system.strands[0].sequence = "TATATA"

    assert len(system.strands) == 3
    assert len(system.nucleotides) == 18
    assert system.strands[0].sequence == "TATATA"
    assert system.strands[2].sequence == "AGAGAG"

    strand_3 = system.strands[1].copy()
    strand_3.sequence = "CCCGGG"
    strand_4 = strand_3.copy()
    strand_4.sequence = "AAATTT"
    system.add_strands({
        0 : strand_3, 
        1 : strand_4
    })

    assert len(system.strands) == 5
    assert len(system.nucleotides) == 30
    print(system.strands)
    assert system.strands[0].sequence == "CCCGGG"
    assert system.strands[1].sequence == "AAATTT"
    assert system.strands[2].sequence == "TATATA"

    print(system)
    print(system.dataframe)
    # system.write_oxDNA()

    system.translate(np.array([10., 0., 0.]))
    

    system.rotate([0., 0., 0., 1.])
    system.rotate([0., 0., 0.])
    system.rotate(
        np.array([
            [1., 0., 0.],
            [0., 1., 0.],
            [0., 0., 0.],
        ])
    )

    system.transform(np.array([
        [1., 0., 0., 0.],
        [0., 1., 0., 0.],
        [0., 0., 1., 0.],
    ]))
    return


if __name__ == "__main__":
    test_System()
